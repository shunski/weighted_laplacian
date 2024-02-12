use std::time::Instant;

use laplacian::obstacles::Obstacle;
use plotters::prelude::*;
use laplacian::WeightedLaplacianHandle;
use laplacian::embedding_lib::Collection;
use alg::lin_alg::ConstVector;

const FOLLOWS_POSITIVE_OF_GRADIENT: bool = true;
const N_HOLES_TO_CREATE: usize = 50;
const STEP_SIZE: f64 = 70000.0;
const DRAW_EDGES: bool = true;

// if COLORING_BY_EIGENVEC is false, then it colors the graph based on the weights
const COLORING_BY_EIGENVEC: bool = false;

fn main() -> Result<(), Box<dyn std::error::Error>>  {
    // let embedding_type = Collection::FourPoints;
    let embedding_type = Collection::Meshlike;
    let mut points = embedding_type.get();
    let original_points = points.clone();
    let n_robots = points.len();
    
    let gradient_type_desc = if FOLLOWS_POSITIVE_OF_GRADIENT { "CollapsingHoles" } else { "CreatingHoles" };
    let file_name = format!("program_outputs/WITH_OBSTACLES==={}==={}===STEP_SIZE:{STEP_SIZE}===EDGES:{:?}===COLORING_BY:{:?}==={}_ROBOTS.jpg", 
        embedding_type.name(), 
        gradient_type_desc, 
        DRAW_EDGES,
        if COLORING_BY_EIGENVEC {"EIGENVEC"} else {"WEIGHTS"},
        points.len()
    );
    let root = BitMapBackend::new(
        &file_name, 
        (1000, 1000)
    ).into_drawing_area();
    root.fill(&WHITE)?;
    let panels = root.split_evenly((3,3));

    let mut obstacles = Obstacle::get();
    for obstacle in &obstacles {
        obstacle.remove_points_contained(&mut points);
    }

    {
        let smallest_d = 4.0;

        for obstacle in &mut obstacles {
            // subdivide the edges so that the delaunay triangulation "captures" the contour of the obstalces.
            obstacle.subdivide_contour(smallest_d);
            // add some noise to the contour so that the points of the contour are not in the general position.
            obstacle.add_noise(0.001);
        }
    };

    for t in 0..=80 {
        let start_time = Instant::now();
        let n_holes = if FOLLOWS_POSITIVE_OF_GRADIENT {None} else {Some(N_HOLES_TO_CREATE)};
        let obstacles_op = if t==0 {Some(obstacles.clone())} else {None};
        let handle = WeightedLaplacianHandle::new_from_points(points, n_holes, obstacles_op);
        let elapsed_time = start_time.elapsed();
        println!("step {t} took {:.5} seconds", elapsed_time.as_secs_f64());

        if t % 10 == 0 { 
            // ------------------ Drawing -------------------------
            let panel = &panels[t / 10];
            let mut graphic = ChartBuilder::on(&panel)
                .margin(10)
                .caption(format!("t={t}"), ("san-serif", 20))
                .build_cartesian_2d(-50.0..50.0, -50.0..50.0)?;
            
            // draw edges
            if DRAW_EDGES {
                
                if COLORING_BY_EIGENVEC {
                    let norm = handle.fiedler_vec_iter().map(|x| x.abs()).max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
                    for (ends, value) in handle.edge_pos_iter().zip(handle.fiedler_vec_iter()) {
                        graphic.draw_series(LineSeries::new(
                            ends.into_iter(),
                            get_color_bwr(value / norm).filled().stroke_width(2)
                        ))?;
                    }
                } else {
                    // let norm = handle.weight_iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap() / 2.0;
                    let mut edge_weights = handle
                        .edge_pos_iter()
                        .zip(handle.weight_iter())
                        .filter(|([x,y],_)| 
                            original_points.contains( &ConstVector::from([x.0,  x.1]) ) && 
                            original_points.contains( &ConstVector::from([y.0,  y.1]) ) 
                        )
                        .collect::<Vec<_>>();
                    edge_weights.sort_by(|(_, w1), (_, w2)| w1.partial_cmp(w2).unwrap() );
                    for (ends, w) in edge_weights {
                        graphic.draw_series(LineSeries::new(
                            ends.into_iter(),
                            get_color_bwr( (w / 2.0).min(1.0) ).filled().stroke_width(2)
                        ))?;
                    }
                }
            }
    
            
            // draw vertices
            graphic.draw_series( handle.vertex_pos_iter(Some(&obstacles))
                .map(|(x,y)| Circle::new((x, y), 2, BLACK.filled()))
            )?;

            // draw obstacles 
            for obstacle in &obstacles {
                graphic.draw_series(std::iter::once(Polygon::new(
                    obstacle.contour_iter().map(|v| (v[0], v[1])).collect::<Vec<_>>(),
                    &YELLOW.mix(0.4)
                )))?;

                graphic.draw_series(std::iter::once(PathElement::new(
                    obstacle.contour_iter()
                        .chain([obstacle.contour_iter().next().unwrap()])
                        .map(|v| (v[0], v[1]))
                        .collect::<Vec<_>>(),
                    &BLACK
                )))?;
            }
        }
        
        
        // move the points by gradient discent
        points = handle.move_points::<FOLLOWS_POSITIVE_OF_GRADIENT>(STEP_SIZE, 2.0, Some(&obstacles));
        println!("");
        root.present()?;
    }

    Ok(())
}

fn get_color_bwr(val: f64) -> RGBColor {
    assert!( -1.0 <= val && val <= 1.0 );
    if val >= 0.0 {
        RGBColor(255, (255.0 * (1.0 - val)) as u8, (255.0 * (1.0 - val)) as u8)
    } else {
        RGBColor((255.0 * (1.0 + val)) as u8, (255.0 * (1.0 + val)) as u8, 255)
    }
}