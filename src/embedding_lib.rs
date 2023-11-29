use std::f64::consts::PI;

use alg::lin_alg::ConstVector;
use geom::delaunay_triangulation;
use rand::prelude::*;

pub enum Collection {
    FourPoints,
    OneHole(usize),
    Meshlike,
    Uniform(usize),
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
                    let n = ((v-w).two_norm() / 2.5497) as usize;
                    let mut i = 0;
                    while i < n {
                        let s = i as f64 / n as f64;
                        let p = v + (w-v)*s + ConstVector::from([rng.gen_range(-2.0..2.0), rng.gen_range(-2.0..2.0)]);

                        if out.iter().all(|&q| (p - q).two_norm() > 1.0) {
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

                let n_random_points = 30;
                let mut i = 0;
                while i < n_random_points {
                    let p = ConstVector::from([rng.gen_range(-35.0..35.0), rng.gen_range(-35.0..35.0)]);

                    if out.iter().all(|&q| (q - p).two_norm() > 1.0 ) {
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
            }
        }
    }

    pub fn name(&self) -> String {
        match self {
            Self::FourPoints => String::from( "FourPoints" ),
            Self::OneHole(n) => format!( "OneHoleBy({n})Points" ),
            Self::Meshlike => format!( "Meshlike" ),
            Self::Uniform(n) => format!("Uniform({n})")
        }
    }
}