use alg::lin_alg::{Matrix, SubMatrix, ConstVector, SparseMatrix};
use geom::delaunay_triangulation;
use analysis::{MultiVarFn, AtomicSmoothFn, Function};
use obstacles::Obstacle;
// use std::time::Instant;

#[allow(unused)]
pub struct WeightedLaplacianHandle {
    points: Vec<ConstVector<f64, 2>>,
    points_vector_form: Matrix<f64>,
    edges: Vec<[usize; 2]>,
    dlambda_dvertices: Matrix<f64>,
    weights_evaluated: Matrix<f64>,
    fiedler_vec: Matrix<f64>,
    laplacian_evaluated: Matrix<f64>,
}

impl WeightedLaplacianHandle {
    pub fn new_from_points(points: Vec<ConstVector<f64, 2>>, n_holes: Option<usize>, obstacles: Option<Vec<Obstacle>>) -> Self {
        let n_holes = n_holes.unwrap_or(1);

        println!("There are {} vertices.", points.len());

        // let start_time = Instant::now();

        // get the delaunay triangulation
        let (_, triangles) = delaunay_triangulation( &points );
        // println!("triangles = {triangles:?}");

        let mut n_points = points.len();

        // then convert points into the form of 2*n dimensional vector 
        let points_vector_form = {
            let mut out = Matrix::zero(points.len() * 2, 1);
            for (i, &x) in points.iter().map(|x| x.iter() ).flatten().enumerate() {
                out[(i,0)] = x;
            }
            out
        };


        // remove triangles that are quotiented out.
        // If 'obstacles' is 'None', then there is no change in 'triangles'.
        // 'removed_triangles' contains removed triangles.
        let (triangles, removed_triangles) = if let Some(obstacles) = obstacles {
            let mut removed_triangles = Vec::new();
            let triangles = triangles.into_iter()
                .map(|t_idx| (t_idx, [points[t_idx[0]], points[t_idx[1]], points[t_idx[2]]]))
                .filter_map(|(t_idx, t)|
                    if obstacles.iter().any(|obstacle| obstacle.intersects(t)) {
                        removed_triangles.push(t_idx);
                        None 
                    } else {
                        Some(t_idx)
                    }
                )
                .collect();
            (triangles, removed_triangles)
        } else {
            (triangles, Vec::new())
        };

        // each element in 'edges' corresponds to two indeces of the boundaries
        let (edges, removed_edges) = {
            let mut edges: Vec<_> = triangles.iter().map(|&[i,j,k]| {
                    [[i,j], [i,k], [j,k]].into_iter()
                })
                .flatten()
                .collect();
            edges.sort();
            edges.dedup();
            

            let mut removed_edges: Vec<_> = removed_triangles.iter().map(|&[i,j,k]| {
                    [[i,j], [i,k], [j,k]].into_iter()
                })
                .flatten()
                .collect();
            removed_edges.sort();
            removed_edges.dedup();

            edges = edges.into_iter().filter(|edge| removed_edges.binary_search(&edge).is_err()).collect();
            (edges, removed_edges)
        };

        // update 'n_points'
        n_points = (0..points.len()).filter(|&i| removed_edges.iter().all(|&e| e[1]!=i ) ).count();
        
        // println!("edges = {edges:?}");
        println!("There are {} edges.", edges.len());

        let n = edges.len();

        // compute 'signs'
        let mut signs = Matrix::zero(n, n);
        for (e1, e2, sign) in triangles.iter()
            .map(|&[i,j,k]|{
                [([i,j], [i,k], -1), ([i,j], [j,k], 1), ([i,k], [j,k],-1)]
            })
            .flatten()
        {
            let i = if let Ok(i) = edges.binary_search(&e1) {
                i
            } else {
                continue;
            };

            let j = if let Ok(j) = edges.binary_search(&e2) {
                j
            } else {
                continue;
            };

            signs[(i,j)] = sign;
            signs[(j,i)] = sign;
        }

        let mut tmp = Matrix::<AtomicSmoothFn>::zero(n, n);

        // compute laplacian
        let laplacian = MultiVarFn::new((n*n, n)).set(|k, f, x| {
            let (i,j) = (k/n, k%n);

            if i > j { *f = tmp[(j,i)].clone(); return };
            
            // now i <= j
            *f = if i==j {
                let mut f = AtomicSmoothFn::Zero;
                for t in triangles.iter()
                    .filter(|&&t| t.contains( &edges[i][0] ) && t.contains( &edges[i][1] ))
                {
                    let v_not_in_edge = *t.iter().find( |v| !edges[i].contains(v) ).unwrap();

                    let mut e1 = [v_not_in_edge, edges[i][0]];
                    e1.sort();
                    let e1 = edges.binary_search(&e1).ok();

                    let mut e2 = [v_not_in_edge, edges[i][1]];
                    e2.sort();
                    let e2 = edges.binary_search(&e2).ok();

                    f += match (e1, e2) {
                        (Some(e1), Some(e2)) => x[e1] * x[e2],
                        (Some(e1), None) => x[e1] * 1.0,
                        (None, Some(e2)) => x[e2] * 1.0,
                        (None, None) => AtomicSmoothFn::Zero,
                    };
                }

                tmp[(i,j)] = f;
                tmp[(i,j)].clone()
            } else if let Some(&common_vertex) = edges[i].iter().find(|k| edges[j].contains(k) ) {
                // if the edges i and j are adjacent, ...
                let mut the_other_edge = [
                    *edges[i].iter().find(|&&k| k != common_vertex ).unwrap(),
                    *edges[j].iter().find(|&&k| k != common_vertex ).unwrap()
                ];
                the_other_edge.sort();

                tmp[(i,j)] = if let Ok(the_other_edge) = edges.binary_search(&the_other_edge) {
                    x[i].powf(0.5) * x[j].powf(0.5) * signs[(i,j)] as f64 * x[the_other_edge]
                } else if removed_edges.binary_search(&the_other_edge).is_ok() {
                    x[i].powf(0.5) * x[j].powf(0.5) * signs[(i,j)] as f64
                } else {
                    AtomicSmoothFn::Zero
                };

                tmp[(i,j)].clone()
            } else {
                AtomicSmoothFn::Zero
            };
        });
        // println!( "laplacian={laplacian}" );
        let dl_dweights = laplacian.clone().jacobian();
        // println!("dl_dweights = {:?}", dl_dweights.size());
        // println!("dl_dweights = {dl_dweights}");

        let (weights, dweights_dvertices) = Self::compute_weights(&edges, points.len());
        // println!("weights = {}", weights);
        // println!("dweights_dvertices = {}", dweights_dvertices);

        let weights_evaluated = weights.eval(&points_vector_form);
        let laplacian_evaluated = {
            let mut out = Matrix::zero( edges.len(), edges.len() );
            for (i, j) in (0..edges.len()).map(|i| (0..edges.len()).map(move |j| (i,j) )).flatten() {
                out[(i,j)] = laplacian[i*edges.len()+j].eval(&weights_evaluated);
            }
            out
        };
        // println!("laplacian_evaluated = {:10}", laplacian_evaluated);
        assert!((0..edges.len()).map(|i| (i+1..edges.len()).map(move |j| (i, j)) ).flatten().all(|(i,j)| laplacian_evaluated[(i,j)] == laplacian_evaluated[(j,i)] ), "laplacian not symmetric");

        // let elapsed_time = start_time.elapsed();
        // println!("constructing laplacian took {:.5} seconds", elapsed_time.as_secs_f64());
        // let start_time = Instant::now();
        let (spectrum, v) = laplacian_evaluated.clone().spectrum_with_n_th_eigenvec_symmetric(n_points-1 + n_holes-1);
        println!("spectrum={}", spectrum);
        // println!("v={v:?}");
        // let elapsed_time = start_time.elapsed();
        // println!("computation of spectrum took {:.5} seconds", elapsed_time.as_secs_f64());

        
        // let start_time = Instant::now();
        let dl_dweights_evaluated = dl_dweights.eval(&weights_evaluated);
        
        let dlambda_dweights = {
            let mut out: Matrix<f64> = Matrix::zero(1, edges.len());
            for col in 0..edges.len() {
                for (row, &x) in dl_dweights_evaluated.col_iter(col) {
                    let (i, j) = (row / edges.len(), row % edges.len());
                    out[(0, col)] += x * v[(i,0)] * v[(j,0)];
                };
            }
            out
        };

        // let dlambda_dweights_ = {
        //     let mut dlambda_dl = Matrix::zero(1, edges.len()*edges.len());
        //     for (k, (i,j)) in (0..edges.len()).map(|i| (0..edges.len()).map(move |j| (i,j) )).flatten().enumerate() {
        //         dlambda_dl[(0,k)] = v[(i,0)] * v[(j,0)];
        //     }
        //     dlambda_dl * dl_dweights_evaluated
        // };
        // assert!((&*dlambda_dweights-dlambda_dweights_).frobenius_norm() / dlambda_dweights.frobenius_norm() < 0.001);
        // println!("dlambda_dweights = {dlambda_dweights:.5}");
        
        let dlambda_dvertices = dlambda_dweights * dweights_dvertices.eval(&points_vector_form);
        // println!("dlambda_dvertices = {:4}", dlambda_dvertices);
        // println!("dlambda_dvertices.two_norm() = {:?}", (&*dlambda_dvertices).transpose().two_norm());

        // let elapsed_time = start_time.elapsed();
        // println!("computation of gradient took {:.5} seconds", elapsed_time.as_secs_f64());

        Self {
            points,
            points_vector_form,
            edges,
            dlambda_dvertices,
            weights_evaluated,
            fiedler_vec: v[(.., 0)].as_matrix(),
            laplacian_evaluated
        }
    }

    fn compute_weights(edges: &[[usize;2]], n_vertices: usize) -> (MultiVarFn, SparseMatrix<AtomicSmoothFn>) {
        let weights = MultiVarFn::new((edges.len(), 2*n_vertices)).set(|i, f, x| {
            let [k, l] = edges[i];
            *f = ((x[2*k]-x[2*l]).powi(2) + (x[2*k+1]-x[2*l+1]).powi(2)).powf(-0.5);
        });
        (weights.clone(), weights.jacobian())
    }

    pub fn edge_pos_iter(&'_ self) -> impl '_ + Iterator<Item = [(f64, f64); 2]> {
        self.edges.iter()
            .map(|&[x, y]| [(self.points[x][0], self.points[x][1]), (self.points[y][0], self.points[y][1])] )
    }

    pub fn weight_iter(&'_ self) -> impl '_ + Iterator<Item = f64> {
        (0..self.edges.len()).map( |i| self.weights_evaluated[(i,0)] )
    }

    pub fn fiedler_vec_iter(&'_ self) -> impl '_ + Iterator<Item = f64> {
        (0..self.edges.len()).map( |i| self.fiedler_vec[(i,0)] )
    }

    pub fn fiedler_vec(&'_ self) -> &'_ SubMatrix<f64> {
        &*self.fiedler_vec
    }

    pub fn vertex_pos_iter(&self) -> impl '_ + Iterator<Item = (f64, f64)> {
        self.points.iter().map(|x| (x[0], x[1]) )
    }

    pub fn move_points<const POSITIVE_GRAD: bool>(self, step_size: f64, max_move: f64) -> Vec<ConstVector<f64, 2>> {
        let sign = if POSITIVE_GRAD { 1.0 } else { -1.0 };
        let dv = {
            let mut dv = self.dlambda_dvertices.transpose() * step_size * sign;
            if dv.two_norm() > max_move {
                dv *= max_move / dv.two_norm();
            }
            dv
        };
        // let mut dv = self.dlambda_dvertices.transpose() * sign;
        // dv /= dv.two_norm(); 
        let new_pos = self.points_vector_form + dv;
        let points = (0..(new_pos.size().0)/2)
            .map(|i| ConstVector::from( [new_pos[(i*2, 0)], new_pos[(i*2+1, 0)]] ))
            .collect::<Vec<_>>();
        
        points
    }
}

pub mod embedding_lib;
pub mod obstacles;

#[cfg(test)]
mod test {
    use alg::matrix;
    use alg::lin_alg::Matrix;
    use alg::lin_alg::ConstVector;
    use geom::delaunay_triangulation;

    use crate::WeightedLaplacianHandle;
    #[test]
    fn from_points() {
        let rt3 = 3_f64.sqrt();
        let rt3_inv_rt = (1.0/rt3).sqrt();
        let b2t = matrix!(f64;
            [[rt3_inv_rt, -rt3_inv_rt,  0.0, 1.0, 0.0, 0.0],
             [rt3_inv_rt,  0.0, -rt3_inv_rt, 0.0, 1.0, 0.0],
             [0.0,  rt3_inv_rt, -rt3_inv_rt, 0.0, 0.0, 1.0]]
        );
        let l = b2t.clone().transpose() * b2t;


        let points = vec![
            ConstVector::from([0.0, 0.0]),
            ConstVector::from([0.0, 1.0]),
            ConstVector::from([-rt3/2.0, -1.0/2.0]),
            ConstVector::from([ rt3/2.0, -1.0/2.0])
        ];
        let handle = WeightedLaplacianHandle::new_from_points(points, None, None);

        assert!(
            ((l.clone() - handle.laplacian_evaluated.clone()).frobenius_norm()/(l.clone()).frobenius_norm()).abs() < 0.000001,
            "l={:.5}, computed l={:.5}", l, handle.laplacian_evaluated
        );
    }

    #[test]
    fn gradient() {
        let rt3 = 3_f64.sqrt();
        let step_size = 0.001;
        let points = vec![
            ConstVector::from([0.0, 0.0]),
            ConstVector::from([0.0, 1.0]),
            ConstVector::from([-rt3/2.0, -1.0/2.0]),
            ConstVector::from([ rt3/2.0, -1.0/2.0])
        ];
        let n = points.len();

        
        // the first computation
        let triangulation1 = delaunay_triangulation(&points).1;
        let handle = WeightedLaplacianHandle::new_from_points(points, None, None);
        let lambda1 = handle.laplacian_evaluated.as_matrix().spectrum_with_n_th_eigenvec_symmetric(n-1).0;
        println!("{:.5}", handle.laplacian_evaluated.as_matrix().spectrum());
        let d_norm = (&*handle.dlambda_dvertices).transpose().two_norm();
        let points = handle.move_points::<true>(step_size, 1.0);

        // the second computation
        let triangulation2 = delaunay_triangulation(&points).1;
        let handle = WeightedLaplacianHandle::new_from_points(points, None, None);
        let lambda2 = handle.laplacian_evaluated.as_matrix().spectrum_with_n_th_eigenvec_symmetric(n-1).0;

        // make sure that the triangulation is the same (so that there is no discrete chnage in the complex)
        assert_eq!(triangulation1, triangulation2);

        assert!(
            ((lambda2 - lambda1 - d_norm*d_norm*step_size)/(d_norm*d_norm*step_size)).abs() < 0.01,
            "lambda2 - lambda1 = {}, but\n d_norm*d_norm*step_size = {}", 
            lambda2 - lambda1, d_norm*d_norm*step_size
        );
    }
}

pub struct WeightedLaplacianHandleCustomWeight {
    weights: Matrix<f64>,
    laplacian: MultiVarFn,
    n_points: usize,
    n_edges: usize,
    triangles_with_low_w: Vec<usize>,
    edges_with_low_w: Vec<usize>,
}

impl WeightedLaplacianHandleCustomWeight {
    pub fn new(points: Vec<ConstVector<f64, 2>>) -> (Self, Vec<[ConstVector<f64, 2>;2]>, Vec<[ConstVector<f64, 2>;3]>) {
        // get the delaunay triangulation
        let (_, mut triangles) = delaunay_triangulation( &points );
        debug_assert!(triangles.iter().all(|[a,b,c]| a < b && b < c ));
        triangles.sort();
        debug_assert!(triangles.iter().zip(triangles.iter().skip(1)).all(|(x,y)| x < y ) );

        // each element in 'edges' corresponds to two indeces of the boundaries
        let edges: Vec<_> = {
            let mut edges: Vec<_> = triangles.iter().map(|&[i,j,k]| {
                    [[i,j], [i,k], [j,k]].into_iter()
                })
                .flatten()
                .collect();

            edges.sort();
            edges.dedup();
            edges
        };

        let n = edges.len();
        let m = triangles.len();

        let mut signs = Matrix::zero(edges.len(), edges.len());
        for (e1, e2, sign) in triangles.iter()
            .map(|&[i,j,k]|{
                [([i,j], [i,k], -1), ([i,j], [j,k], 1), ([i,k], [j,k],-1)].into_iter()
            })
            .flatten()
        {
            let i = edges.binary_search(&e1).unwrap();
            let j = edges.binary_search(&e2).unwrap();
            signs[(i,j)] = sign;
            signs[(j,i)] = sign;
        }

        let mut tmp = Matrix::<AtomicSmoothFn>::zero(n, n);

        // compute laplacian
        let laplacian = MultiVarFn::new((n*n, n+m)).set(|k, f, x| {
            let (i,j) = (k/n, k%n);

            if i > j { *f = tmp[(j,i)].clone(); return };
            
            // now i <= j
            *f = if i==j {
                let the_other_two_edges = triangles.iter()
                    .filter(|&&t| t.contains( &edges[i][0] ) && t.contains( &edges[i][1] ))
                    .map(|t| {
                        let v_not_in_edge = *t.iter().find( |v| !edges[i].contains(v) ).unwrap();

                        let mut e1 = [v_not_in_edge, edges[i][0]];
                        e1.sort();
                        let e1 = edges.binary_search(&e1).unwrap();

                        let mut e2 = [v_not_in_edge, edges[i][1]];
                        e2.sort();
                        let e2 = edges.binary_search(&e2).unwrap();

                        let mut face = [v_not_in_edge, edges[i][0], edges[i][1]];
                        face.sort();
                        let face = triangles.binary_search(&face).unwrap();

                        [e1, e2, face]
                    } )
                    .take(2)
                    .collect::<Vec<_>>();



                // each edge must have either one or two incident face
                assert!( the_other_two_edges.len() == 1 || the_other_two_edges.len() == 2 );
                let mut g = AtomicSmoothFn::Zero;
                for [a,b, face] in the_other_two_edges {
                    g += x[n+face].powi(2) / x[i].powi(2);
                } 

                tmp[(i,j)] = g;
                tmp[(i,j)].clone()
            } else if edges[i].iter().find(|k| edges[j].contains(k) ).is_some() {
                // if the edges i and j are adjacent, ...
                let face = {
                    let mut face = vec![ edges[i][0], edges[i][1], edges[j][0], edges[j][1]];
                    face.sort();
                    face.dedup();
                    let face: [usize; 3] = face.try_into().unwrap();
                    triangles.binary_search(&face)
                };

                tmp[(i,j)] = if let Ok(face) = face {
                    x[n+face].powi(2) / x[i] / x[j] * signs[(i,j)] as f64
                } else {
                    AtomicSmoothFn::Zero
                };

                tmp[(i,j)].clone()
            } else {
                AtomicSmoothFn::Zero
            };
        });

        let mut weights = Matrix::zero(n+m, 1);
        for i in 0..n+m {
            weights[(i, 0)] = 1.0;
        }

        let triangles_with_low_w = (0..triangles.len()).filter(|&i| triangles[i].contains(&15) || triangles[i].contains(&45) ).collect::<Vec<_>>();
        let edges_with_low_w = (0..edges.len()).filter(|&i| edges[i].contains(&15) || edges[i].contains(&45) ).collect::<Vec<_>>();

        let out = Self {
            weights,
            laplacian,
            n_points: points.len(),
            n_edges: edges.len(),
            triangles_with_low_w,
            edges_with_low_w
        };

        let edges: Vec<_> = edges.into_iter().map(|[i,j]| [points[i], points[j]] ).collect();
        let triangles: Vec<_> = triangles.into_iter().map(|[i,j, k]| [points[i], points[j], points[k]] ).collect();

        (out, edges, triangles)
    }
}

impl Iterator for WeightedLaplacianHandleCustomWeight {
    type Item = (Matrix<f64>, Matrix<f64>, Matrix<f64>);

    fn next(&mut self) -> Option<Self::Item> {
        for idx in &self.triangles_with_low_w {
            self.weights[(self.n_edges + idx, 0)] /= 1.5;
        }

        for idx in &self.edges_with_low_w {
            self.weights[(*idx, 0)] /= 1.5;
        }

        
        let l = self.laplacian.eval(&self.weights);
        let mut laplacian = Matrix::zero(self.n_edges, self.n_edges);
        for (i,j) in (0..self.n_edges).map(|i| (0..self.n_edges).map(move |j| (i,j) ) ).flatten() {
            laplacian[(i,j)] = l[(i*self.n_edges + j, 0)];
        }
        assert!((0..self.n_edges).all(|i| (i..self.n_edges).all( |j| laplacian[(i,j)] == laplacian[(j,i)] )));
        let (spectrum, v) = laplacian.spectrum_with_n_smallest_eigenvecs_symmetric(self.n_points+10-1);
        let spectrum = spectrum[(self.n_points-2.., 0)].as_matrix();
        let v = v[(..,self.n_points-2..)].as_matrix();

        println!("n_points = {}", self.n_points);
        println!("self.n_edges = {}", self.n_edges);
        // println!("spectrum={spectrum:?}");

        Some((self.weights.clone(), spectrum, v)) 
    }
}