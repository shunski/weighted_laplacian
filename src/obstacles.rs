use std::f64::consts::PI;

use alg::lin_alg::ConstVector;

#[derive(Clone)]
pub struct Obstacle {
    center: ConstVector<f64, 2>,
    contour: Vec<ConstVector<f64, 2>>,
}


impl Obstacle {
    pub fn get() -> Vec<Self> {
        let mut out = Vec::new();
    
        let resolution = 100;
        let diameter = 10.0;
        let center = ConstVector::from([0.0, 0.0]);
    
        let contour = (0..resolution).map( |i| {
            let t = 2.0 * PI / resolution as f64 * i as f64;
            ConstVector::from([t.cos(), t.sin()]) * diameter + center
        }).collect::<Vec<_>>();
        let obstacle = Obstacle{center, contour};
    
        out.push(obstacle);
    
        out
    }

    pub fn contour_iter(&'_ self) -> impl Iterator<Item = &'_ ConstVector<f64, 2>> {
        self.contour.iter().skip(1).chain([self.contour.first().unwrap()])
    }
    
    
    pub fn contains(&self, p: &ConstVector<f64, 2>) -> bool {
        let p = *p;
        self.contour.iter()
            .zip( self.contour.iter().skip(1).chain([self.contour.first().unwrap()]) )
            .all( |(&v, &w)| {
                let v_w = w - v;
                let perp = ConstVector::from([-v_w[1], v_w[0]]);
                debug_assert_eq!(v_w.dot(perp), 0.0, "computation of 'perp' is wrong.");
                (self.center-v).dot(perp) * (p-v).dot(perp) > 0_f64
            } )
    }

    pub fn intersects(&self, triangle: [ConstVector<f64, 2>; 3]) -> bool {
        let edges = [
            (triangle[0], [triangle[1], triangle[2]]),
            (triangle[1], [triangle[0], triangle[2]]),
            (triangle[2], [triangle[0], triangle[1]])
        ];
    
        self.contour.iter().any(|&p| {
            edges.into_iter().all( |(v, e)| {
                let e_perp = ConstVector::from([-(e[1]-e[0])[1], (e[1]-e[0])[0]]);
                let e0_v = v-e[0];
                let e0_p = p-e[0];
                e0_v.dot(e_perp) * e0_p.dot(e_perp) > 0_f64
            })
        })
    }

    pub fn remove_points_contained(&self, points: &mut Vec<ConstVector<f64, 2>>) {
        *points = points.iter()
            .copied()
            .filter(|&p| !self.contains(&p) )
            .collect::<Vec<_>>();
    }
}