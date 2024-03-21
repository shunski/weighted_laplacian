use plotters::prelude::*;
use laplacian::WeightedLaplacianHandle;
use laplacian::embedding_lib::Collection;

fn main() -> Result<(), Box<dyn std::error::Error>>  {
    let embedding_type = Collection::ThreeHoles;
    let points = embedding_type.get();

    let handle = WeightedLaplacianHandle::new_from_points(points, None, None);
    
    
    let weight_file_name = format!("program_outputs/StaticNetwork===Weight==={}.jpg", 
        embedding_type.name()
    );

    let histogram_file_name = format!("program_outputs/StaticNetwork===Weight==={}.jpg", 
        embedding_type.name()
    );

    let eigenvec_file_name = format!("program_outputs/StaticNetwork===Weight==={}.jpg", 
        embedding_type.name()
    );

    let weight_root = BitMapBackend::new(
        &weight_file_name, 
        (500, 500)
    ).into_drawing_area();
    weight_root.fill(&WHITE)?;

    let histogram_root = BitMapBackend::new(
        &histogram_file_name, 
        (500, 500)
    ).into_drawing_area();
    histogram_root.fill(&WHITE)?;

    let eigenvec_root = BitMapBackend::new(
        &eigenvec_file_name, 
        (1500, 500)
    ).into_drawing_area();
    weight_root.fill(&WHITE)?;
    let eigenvec_panels = eigenvec_root.split_evenly((1,3));



    // ======================== DRAWING HISTOGRAM =========================================
    let mut histogram = ChartBuilder::on(&histogram_root)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .margin(5)
        .build_cartesian_2d((1u32..10u32).into_segmented(), 0f64..1.0)?;

    histogram
        .configure_mesh()
        .disable_mesh()
        .axis_desc_style(("sans-serif", 10))
        .draw()?;

    histogram.draw_series(
        Histogram::vertical(&histogram)
            .style(BLACK)
            .data(handle.spectrum_iter() ),
    )?;


    // ========================= DRAWING THE WEIGHTS =================================
    let mut w = ChartBuilder::on(&weight_root)
        .margin(10)
        .build_cartesian_2d(-50.0..50.0, -50.0..50.0)?;

    // draw triangles
    let weight_iter = (0..weights.size().0).map(|i| weights[(i+edges.len(), 0)] );
    for (t, weight) in handle.triangles.iter().zip( weight_iter ) {
        w.draw_series(std::iter::once(Polygon::new(
            t.iter().map(|x| (x[0], x[1]) ).collect::<Vec<_>>(),
            get_color_bwr( weight )
        )))?;
    }

    // draw edges
    let weight_iter = (0..weights.size().0).map(|i| weights[(i,0)] );
    for (ends, weight) in edges.iter().zip( weight_iter ) {
        w.draw_series(LineSeries::new(
            ends.into_iter().map(|v| (v[0], v[1]) ),
            get_color_bwr( weight ).filled().stroke_width(2)
        ))?;
    }

    // draw vertices
    w.draw_series( 
        points.iter().map(|v| Circle::new((v[0], v[1]), 2, BLACK.filled()))
    )?;



    // ========================= DRAWING THE SECOND EIGENVECTOR =================================
    let mut v2 = ChartBuilder::on(&v2)
    .margin(10)
    .build_cartesian_2d(-50.0..50.0, -50.0..50.0)?;

    // draw edges
    let norm = (0..edges.len()).map(|i| eigenvec[(i,1)].abs() ).max_by(|x,y| x.partial_cmp(y).unwrap() ).unwrap() + 0.001;
    for (i, edge) in edges.iter().enumerate() {
        v2.draw_series(LineSeries::new(
            edge.into_iter().map(|v| (v[0], v[1]) ),
            get_color_bwr( eigenvec[(i,1)]/norm ).filled().stroke_width(2)
        ))?;
    }

    // draw vertices
    v2.draw_series( 
        points.iter().map(|v| Circle::new((v[0], v[1]), 2, BLACK.filled()))
    )?;
    

    histogram_root.present()?;
    weight_root.present()?;
    eigenvec_root.present()?;

    Ok(())
}

fn get_color_bwr(val: f64) -> RGBColor {
    assert!( -1.0 <= val && val <= 1.0, "val == {val}");
    if val >= 0.0 {
        RGBColor(255, (255.0 * (1.0 - val)) as u8, (255.0 * (1.0 - val)) as u8)
    } else {
        RGBColor((255.0 * (1.0 + val)) as u8, (255.0 * (1.0 + val)) as u8, 255)
    }
}