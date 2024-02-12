use std::f64::consts::PI;
use alg::lin_alg::ConstVector;
use rand::Rng;

#[derive(Clone)]
pub struct Obstacle {
    contour: Vec<ConstVector<f64, 2>>,
}


impl Obstacle {
    pub fn get() -> Vec<Self> {
        let mut out = Vec::new();
    
        let obstacle = {
            let contour = vec![
                ConstVector::from([-30.0, -30.0]),
                ConstVector::from([-30.0,  10.0]),
                ConstVector::from([-25.0,  10.0]),
                ConstVector::from([-25.0, -25.0]),
                ConstVector::from([ 10.0, -25.0]),
                ConstVector::from([ 10.0, -30.0]),
            ];
            Obstacle{ contour }
        };
        out.push( obstacle );



        let obstacle = {
            let contour = vec![
                ConstVector::from([-10.0, -5.0]),
                ConstVector::from([-10.0,  0.0]),
                ConstVector::from([ 30.0,  0.0]),
                ConstVector::from([ 30.0, -5.0]),
            ];
            Obstacle{ contour }
        };
        out.push( obstacle );


        let obstacle = {
            let contour = vec![
                ConstVector::from([-10.0, 20.0]),
                ConstVector::from([-10.0, 25.0]),
                ConstVector::from([ 30.0, 25.0]),
                ConstVector::from([ 30.0, 20.0]),
            ];
            Obstacle{ contour }
        };
        out.push( obstacle );
    
        out
    }
    

    pub fn contour_iter(&'_ self) -> impl Iterator<Item = &'_ ConstVector<f64, 2>> {
        self.contour.iter()
    }

    pub fn contour_size(&'_ self) -> usize {
        self.contour.len()
    }
    
    
    pub fn contains(&self, p: &ConstVector<f64, 2>) -> bool {
        let p = *p;
        let integral = self.contour.iter()
            .zip( self.contour.iter().skip(1).chain([self.contour.first().unwrap()]) )
            .map( |(&v, &w)| {
                let v = v - p;
                let w = w - p;
                let v_perp = ConstVector::from([-v[1], v[0]]);
                assert!(v.dot(v_perp)==0_f64);

                let a = v.two_norm();
                let b = w.two_norm();
                let c = (w-v).two_norm();

                let sign = if (w-v).dot(v_perp)>0.0 {1.0} else {-1.0};

                let theta = sign * ((a*a + b*b - c*c)/(2.0*a*b)).acos();
                theta
            } )
            .sum::<f64>();
        
        // if the integral is zero, then 'p' is not contained in self.
        // if the integral is 2*PI (or -2*PI), then 'p' is contained in self.
        integral.abs() > PI
    }

    pub fn intersects(&self, triangle: [ConstVector<f64, 2>; 3]) -> bool {
        let edges = [
            (triangle[0], triangle[1]),
            (triangle[1], triangle[0]),
            (triangle[2], triangle[0])
        ];
    
        self.contour.iter()
            .zip( self.contour.iter().skip(1).chain([self.contour.first().unwrap()]) )
            .any(|(&p, &q)| {
                edges.into_iter().any( |(v, w)| {
                    let quad = Self{contour: vec![p,v,q,w]};
                    let x = (p+q)/2.0;
                    let y = (v+w)/2.0;
                    quad.contains(&x) && quad.contains(&y)
                })
        })
    }

    pub fn remove_points_contained(&self, points: &mut Vec<ConstVector<f64, 2>>) {
        *points = points.iter()
            .copied()
            .filter(|&p| !self.contains(&p) )
            .collect::<Vec<_>>();
    }


    pub fn subdivide_contour(&mut self, max_length: f64 ) {
        let mut subdivided_contour = Vec::new();
        for (&p, &q) in self.contour.iter().zip(self.contour.iter().skip(1)) {
            let v = q-p;
            let n = (v.two_norm() / max_length) as usize + 1;
            let u = v / n as f64;
            for i in 0..n {
                subdivided_contour.push(u * i as f64 + p);
            }
        }

        subdivided_contour.push(*self.contour.last().unwrap());

        self.contour = subdivided_contour;
    }

    pub fn add_noise(&mut self, amp: f64) {
        let mut rng = rand::thread_rng();
        for p in &mut self.contour {
            let theta = rng.gen_range(0.0..PI);
            let r = rng.gen_range(-amp..amp);
            *p += ConstVector::from([r*theta.cos(), r*theta.sin()]);
        }
    }
}