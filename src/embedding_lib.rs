use std::f64::consts::{PI, TAU};

use alg::lin_alg::ConstVector;
use geom::delaunay_triangulation;
use rand::prelude::*;


pub enum Collection {
    FourPoints,
    OneHole(usize),
    Meshlike,
    Uniform(usize),
    HexagonalLattice,
    ForObstacles,
    ThreeHoles,
}

impl Collection {
    pub fn get(&self) -> Vec<ConstVector<f64, 2>> {
        let mut rng = rand::thread_rng();
        match self {
            Self::FourPoints => vec![
                ConstVector::from([0_f64,  40.0]),
                ConstVector::from([0_f64, -40.0]),
                ConstVector::from([40.0, 0_f64]),
                ConstVector::from([-40.0, 0_f64])
            ],
            Self::OneHole(n) =>
                (0..*n).map( |i| {
                    let theta = (i as f64 / *n as f64) * 2.0 * PI;
                    let unit_length = 100.0;
                    let unit_length_with_noize = ConstVector::from([theta.cos() + rng.gen_range(-0.1..0.1), theta.sin() + rng.gen_range(-0.1..0.1)])* unit_length;
                    if rng.gen_range(0_f64..1.0) < 0.85 {
                        unit_length_with_noize * 0.7
                    } else {
                        let s = rng.gen_range(0_f64..0.9).sqrt();
                        unit_length_with_noize * 0.7 * s
                    }
                }).collect::<Vec<_>>(),

            Self::Meshlike => {
                let mut out = Vec::new();
                let mut add_points = |[v, w]: [ConstVector<f64, 2>; 2]| {
                    let n = ((v-w).two_norm() / 2.5) as usize;
                    let mut i = 0;
                    while i < n {
                        let s = i as f64 / n as f64;
                        let p = v + (w-v)*s + ConstVector::from([rng.gen_range(-2.0..2.0), rng.gen_range(-2.0..2.0)]);

                        if out.iter().all(|&q| (p - q).two_norm() > 1.5) {
                            out.push( p );
                            i+=1;
                        }
                    }
                };

                let points = vec![
                    ConstVector::from([-30.0,  30.0]),
                    ConstVector::from([-30.0, -35.0]),
                    ConstVector::from([-15.0, -15.0]),
                    ConstVector::from([-40.0, -10.0]),
                    ConstVector::from([ 30.0, -35.0]),
                    ConstVector::from([ 35.0,  35.0]),
                    ConstVector::from([ 20.0, -15.0]),
                    ConstVector::from([  0.0,  40.0]),
                    ConstVector::from([  5.0,  15.0]),
                    ConstVector::from([ 40.0,   0.0]),
                    ConstVector::from([-10.0, -40.0])
                ];

                let (_, triangles) = delaunay_triangulation(&points);
                let edges = {
                    let mut edges = triangles.into_iter()
                        .map(|[a,b,c]| [[a,b], [a,c], [b,c]].into_iter() )
                        .flatten()
                        .collect::<Vec<_>>();
                    edges.sort();
                    edges.dedup();
                    edges
                };

                for [v, w] in edges.into_iter().map(|[i, j]| [points[i], points[j]] ) {
                    add_points([v, w]);
                }

                let n_random_points = 35;
                let mut i = 0;
                while i < n_random_points {
                    let p = ConstVector::from([rng.gen_range(-35.0..35.0), rng.gen_range(-35.0..35.0)]);

                    if out.iter().all(|&q| (q - p).two_norm() > 2.0 ) {
                        out.push(p);
                        i+=1;
                    }
                }

                println!("the complex has {} points", out.len());

                out
            },
            Self::Uniform(n) => {
                let mut out = Vec::new();
                while out.len() < *n {
                    let p = ConstVector::from( [ rng.gen_range(-45.0..45.0), rng.gen_range(-45.0..45.0) ]);
                    if out.iter().all(|&q| (p-q).two_norm() > 3.0 ) {
                        out.push(p);
                    }
                }
                out
            },
            Self::HexagonalLattice => {
                let lattice_basis = [
                    ConstVector::from( [1.0, 0f64]) * 12.0,
                    ConstVector::from( [0.5, 3f64.sqrt()/2.0]) * 12.0,
                ];

                (-10..10).map(|x| (-10..10).map(move |y| (x,y) ))
                    .flatten()
                    .map(|(x,y)| lattice_basis[0] * x as f64 + lattice_basis[1] * y as f64)
                    .filter(|v| v.two_norm() <= 50.0 )
                    .collect::<Vec<_>>()
            },
            Self::ForObstacles => {
                let n = 500;
                let mut out = Vec::new();

                // Randomly add 500 points.
                while out.len() < n {
                    let p = ConstVector::from( [ rng.gen_range(-37.0..37.0), rng.gen_range(-37.0..37.0) ]);
                    if out.iter().all(|&q| (p-q).two_norm() > 2.5 ) {
                        out.push(p);
                    }
                }

                // Then add more points on sides.
                let n_additional = 80;
                while out.len() < n + n_additional {
                    let p = ConstVector::from( [ rng.gen_range(-37.0..37.0), rng.gen_range(-37.0..37.0) ]);
                    let forbidden_range = -32.0..32.0;
                    if out.iter().all(|&q| (p-q).two_norm() > 2.0 ) && !(forbidden_range.contains(&p[0]) && forbidden_range.contains(&p[1])) {
                        out.push(p);
                    }
                }
                
                let circles = [
                    (ConstVector::from([23.0, 10.0]), 8.0),
                    (ConstVector::from([2.0, 10.0]), 8.0),
                    (ConstVector::from([12.0, -17.0]), 10.0),
                    (ConstVector::from([-15.0, -15.0]), 10.0),
                    (ConstVector::from([-16.0, 17.0]), 10.0),
                ];

                for (center, r) in circles {
                    for p in &mut out{
                        let d = (*p - center).two_norm();
                        if d >= r {continue;}

                        let s = (d/r).cbrt();
                        *p = center + (*p - center) * (s * r / d)
                    }
                }

                
                // Remove points that are too close to each other.
                let mut i=0;
                while i < out.len() {
                    let mut j = i+1;
                    while j < out.len() {
                        if (out[i]-out[j]).two_norm() < 1.2 {
                            out.remove(j);
                        } else {
                            j += 1;
                        }
                    }
                    i += 1;
                }

                println!("There are {} robots.", out.len());

                out
            },

            Self::ThreeHoles => {
                let circles = [
                    (ConstVector::from([0.0, 18.0]), 18.0),
                    (ConstVector::from([-16.0, -12.0]), 18.0),
                    (ConstVector::from([16.0, -11.0]), 15.0),
                ];

                let mut out = Vec::new();
                
                while out.len() < 180 {
                    let (center, r) = circles[rng.gen_range(0..3)];
                    let angle = rng.gen_range(0.0..TAU);
                    let noise = rng.gen_range(0_f64..1.0).powi(2) * 4.0 * if rng.gen_bool(0.5) {-1.0} else {1.0};
                    let p = ConstVector::from( [ angle.cos(), angle.sin() ]) * (r + noise) + center;
                    if out.iter().all(|q| (p-*q).two_norm() > 1.5) {
                        out.push(p);
                    }
                }

                while out.len() < 190 {
                    let p = ConstVector::from([rng.gen_range(-40.0..40.0), rng.gen_range(-40.0..40.0)] );
                    if out.iter().all(|q| (p-*q).two_norm() > 7.0) && circles.iter().any(|&(c, r)| (p-c).two_norm() < r ) {
                        out.push(p);
                    }
                }

                out.push(ConstVector::from([-24.0,15.0]));
                out.push(ConstVector::from([-28.0,10.0]));
                out.push(ConstVector::from([0.0,-30.0]));
                out.push(ConstVector::from([3.0,-28.0]));
                out.push(ConstVector::from([26.0,13.0]));
                out.push(ConstVector::from([29.0, 8.0]));

                out
            }
        }
    }

    pub fn name(&self) -> String {
        match self {
            Self::FourPoints => String::from( "FourPoints" ),
            Self::OneHole(n) => format!( "OneHoleBy({n})Points" ),
            Self::Meshlike => format!( "Meshlike" ),
            Self::Uniform(n) => format!("Uniform({n})"),
            Self::HexagonalLattice => String::from( "HexagonalLattice" ),
            Self::ForObstacles => String::from( "ForObstacles" ),
            Self::ThreeHoles => String::from("Three holes")
        }
    }
}