use std::time::Instant;

use plotters::prelude::*;
use laplacian::WeightedLaplacianHandle;
use laplacian::embedding_lib::Collection;

const FOLLOWS_POSITIVE_OF_GRADIENT: bool = true;
const STEP_SIZE: f64 = 100000.0;

fn main() -> Result<(), Box<dyn std::error::Error>>  {
    let embedding_type = Collection::Meshlike;
    let mut points = embedding_type.get();

    for t in 1..=100 {
        let start_time = Instant::now();
        let handle = WeightedLaplacianHandle::new_from_points(points, None, None);
        let elapsed_time = start_time.elapsed();
        println!("step {t} took {:.5} seconds", elapsed_time.as_secs_f64());

        // ------------------ Drawing -------------------------
        let gradient_type_desc = if FOLLOWS_POSITIVE_OF_GRADIENT { "CollapsingHoles" } else { "CreatingHoles" };
        let file_name = format!("program_outputs/{}==={}/{t}.jpg", embedding_type.name(), gradient_type_desc);
        let root = BitMapBackend::new(
            &file_name, 
            (1000, 1000)
        ).into_drawing_area();

        root.fill(&WHITE)?;

        let mut graphic = ChartBuilder::on(&root)
            .caption(format!("t={t}"), ("san-serif", 30))
            .build_cartesian_2d(-50.0..50.0, -50.0..50.0)?;
        
        // draw edges
        let norm = handle.fiedler_vec_iter().map(|x| x.abs()).max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
        for (ends, value) in handle.edge_pos_iter().zip(handle.fiedler_vec_iter()) {
            graphic.draw_series(LineSeries::new(
                ends.into_iter(),
                get_color(value / norm).filled().stroke_width(2)
            ))?;
        }

        // draw vertices
        graphic.draw_series( handle.vertex_pos_iter()
            .map(|(x,y)| Circle::new((x, y), 4, BLACK.filled()))
        )?;
        
        root.present()?;


        // move the points by gradient discent
        points = handle.move_points::<FOLLOWS_POSITIVE_OF_GRADIENT>(STEP_SIZE, 2.0);
    }

    Ok(())
}

fn get_color(val: f64) -> RGBColor {
    assert!( -1.0 <= val && val <= 1.0 );
    if val >= 0.0 {
        RGBColor(255, (255.0 * (1.0 - val)) as u8, (255.0 * (1.0 - val)) as u8)
    } else {
        RGBColor((255.0 * (1.0 + val)) as u8, (255.0 * (1.0 + val)) as u8, 255)
    }
}