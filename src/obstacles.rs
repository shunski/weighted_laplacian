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
                ConstVector::from([-5.0, -5.0]),
                ConstVector::from([-5.0,  0.0]),
                ConstVector::from([ 30.0,  0.0]),
                ConstVector::from([ 30.0, -5.0]),
            ];
            Obstacle{ contour }
        };
        out.push( obstacle );


        let obstacle = {
            let contour = vec![
                ConstVector::from([-5.0, 20.0]),
                ConstVector::from([-5.0, 25.0]),
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

    pub fn contour_segment_iter(&'_ self) -> impl Iterator<Item = (&'_ ConstVector<f64, 2>, &'_ ConstVector<f64, 2>)> {
        self.contour_iter().zip(
            self.contour_iter().skip(1).chain([self.contour.first().unwrap()])
        )
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


    pub fn remove_points_closer_than(&self, threshold: f64, points: &mut Vec<ConstVector<f64, 2>>) {
        *points = points.iter()
            .copied()
            .filter(|&p| self.distance_from(p).two_norm() > threshold )
            .collect::<Vec<_>>();
    }


    pub fn subdivide_contour(&mut self, max_length: f64 ) {
        let mut subdivided_contour = Vec::new();
        for (&p, &q) in self.contour.iter().zip(self.contour.iter().skip(1).chain([self.contour.first().unwrap()])) {
            let v = q-p;
            let n = (v.two_norm() / max_length) as usize + 1;
            let u = v / n as f64;
            for i in 0..n {
                subdivided_contour.push(u * i as f64 + p);
            }
        }

        self.contour = subdivided_contour;
    }

    pub fn add_noise(&mut self, amp: f64) {
        let mut rng = rand::thread_rng();
        for p in &mut self.contour {
            let theta = rng.gen_range(0.0..2.0*PI);
            let r = rng.gen_range(-amp..amp);
            *p += ConstVector::from([r*theta.cos(), r*theta.sin()]);
        }
    }

    // returns a distance from the point to the obstacle as a vector pointing TO the obstacle FROM the point.
    // The length of the returned vector is the distance between the obstacle and the point.
    pub fn distance_from(&self, p: ConstVector<f64, 2>) -> ConstVector<f64, 2> {
        self.contour_segment_iter()
            .map(
                |(&x, &y)| {
                    let v = y - x;
                    let w = p - x;

                    let det = (v[0]*w[1]-w[0]*v[1]).abs();

                    // The height of a parallelogram is its area divided by the base.
                    let d = det / v.two_norm();
                    
                    let n = ConstVector::from([-v[1], v[0]]);

                    // detemine the direction of the output so that it points TO the obstacle FROM the point.
                    let sign = if n.dot(w)>0.0 {-1.0} else {1.0}; 
                    n * (sign*d/n.two_norm())
                } 
            )
            .min_by(|x,y| x.two_norm().partial_cmp( &y.two_norm() ).unwrap())
            .unwrap()
    }
}