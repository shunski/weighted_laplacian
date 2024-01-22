use alg::lin_alg::{Matrix, SubMatrix, ConstVector, SparseMatrix};
use geom::delaunay_triangulation;
use analysis::{MultiVarFn, AtomicSmoothFn};
use analysis::Function;
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
    pub fn new_from_points(points: Vec<ConstVector<f64, 2>>, n_holes: Option<usize>, obstacles: Option<Vec<Vec<ConstVector<f64, 2>>>>) -> Self {
        let n_holes = n_holes.unwrap_or(1);

        println!("There are {} vertices.", points.len());

        // let start_time = Instant::now();

        // get the delaunay triangulation
        let (_, triangles) = delaunay_triangulation( &points );
        // println!("triangles = {triangles:?}");

        // then convert points into the form of 2*n dimensional vector 
        let points_vector_form = {
            let mut out = Matrix::zero(points.len() * 2, 1);
            for (i, &x) in points.iter().map(|x| x.iter() ).flatten().enumerate() {
                out[(i,0)] = x;
            }
            out
        };

        // each element in 'edges' corresponds to two indeces of the boundaries
        let edges: Vec<[usize; 2]> = {
            let mut edges: Vec<_> = triangles.iter().map(|&[i,j,k]| {
                    [[i,j], [i,k], [j,k]].into_iter()
                })
                .flatten()
                .collect();

            edges.sort();
            edges.dedup();
            edges
        };
        // println!("edges = {edges:?}");
        println!("There are {} edges.", edges.len());

        let n = edges.len();

        // compute 'signs'
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
        let laplacian = MultiVarFn::new((n*n, n)).set(|k, f, x| {
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

                        [e1, e2]
                    } )
                    .take(2)
                    .collect::<Vec<_>>();

                // each edge must have either one or two incident face
                assert!( the_other_two_edges.len() == 1 || the_other_two_edges.len() == 2 );
                let mut f = AtomicSmoothFn::Zero;
                for [a,b] in the_other_two_edges {
                    f += x[a] * x[b];
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
        let (_spectrum, v) = laplacian_evaluated.clone().spectrum_with_n_th_eigenvec_symmetric(points.len()-1 + n_holes-1);
        // println!("spectrum={}", spectrum[(points.len()-2, 0)]);
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