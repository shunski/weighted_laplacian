use plotters::prelude::*;
use laplacian::embedding_lib::Collection;
use geom;
use alg::lin_alg::{ConstVector, Matrix};

fn main() -> Result<(), Box<dyn std::error::Error>>  {
    let embedding_type = Collection::ThreeHoles;
    let points = embedding_type.get();
    
    
    let weight_file_name = format!("program_outputs/StaticNetwork===WEIGHT==={}.jpg", 
        embedding_type.name()
    );

    let histogram_file_name = format!("program_outputs/StaticNetwork===SPECTRUM==={}.jpg", 
        embedding_type.name()
    );

    let eigenvec_file_name = format!("program_outputs/StaticNetwork===EIGENVEC==={}.jpg", 
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
    eigenvec_root.fill(&WHITE)?;
    let eigenvec_panels = eigenvec_root.split_evenly((1,3));

    // ======================= Computations ==============================================

    let (_, triangles) = geom::delaunay_triangulation(&points);

    let edges = {
        let mut edges = triangles.iter()
            .flat_map(|&[i,j,k]| [[i,j], [i,k], [j,k]] )
            .collect::<Vec<_>>();

        edges.sort();
        edges.dedup();
        edges
    };


    let b1 = {
        let mut b1 = Matrix::zero( points.len(), edges.len() );
        for (e, &[v,w]) in edges.iter().enumerate() {
            b1[(w,e)] = 1.0; b1[(v,e)] = -1.0;
        }
        b1
    };

    let b2 = {
        let mut b2 = Matrix::zero( edges.len(), triangles.len() );
        for (t, &[v,w, x]) in triangles.iter().enumerate() {
            let e1 = edges.binary_search(&[v,w]).unwrap();
            let e2 = edges.binary_search(&[v,x]).unwrap();
            let e3 = edges.binary_search(&[w,x]).unwrap();
            b2[(e1,t)] = 1.0; 
            b2[(e2,t)] = -1.0; 
            b2[(e3,t)] = 1.0;
        }
        b2
    };


    let mut w1 = Matrix::zero( edges.len(), edges.len() );
    let mut w1_inv = Matrix::zero( edges.len(), edges.len() );
    for (i, e) in edges.iter().enumerate() {
        w1[(i,i)] = w(e, &points);
        w1_inv[(i,i)] = 1.0 / w(e, &points);
    }

    let mut w2_sq = Matrix::zero( triangles.len(), triangles.len() );
    for (i, &[s, t, u]) in triangles.iter().enumerate() {
        let e1 = &edges[edges.binary_search(&[s,t]).unwrap()];
        let e2 = &edges[edges.binary_search(&[s,u]).unwrap()];
        let e3 = &edges[edges.binary_search(&[t,u]).unwrap()];
        w2_sq[(i,i)] = (w(e1, &points) * w(e2, &points) * w(e3, &points)).powi(2);
    }

    let b1 = &*b1;
    let b2 = &*b2;
    let w1 = &*w1;
    let w1_inv = &*w1_inv;
    let w2_sq = &*w2_sq;

    let l = w1*b1.transpose()*b1*w1 + w1_inv*b2*w2_sq*b2.transpose()*w1_inv;

    let (spectrum, eigenvecs) = l.spectrum_with_n_smallest_eigenvecs_symmetric(10);

    let spectrum_iter = (0..10).map(|i| (i as u32 + 1, spectrum[(i,0)]) );

    // ======================== DRAWING HISTOGRAM =========================================
    let max = spectrum[(9,0)];
    let mut histogram = ChartBuilder::on(&histogram_root)
        .x_label_area_size(35)
        .margin(5)
        .build_cartesian_2d((1u32..10u32).into_segmented(), 0f64..max*1.1)?;

    histogram
        .configure_mesh()
        .disable_mesh()
        .axis_desc_style(("sans-serif", 10))
        .draw()?;

    histogram.draw_series(
        Histogram::vertical(&histogram)
            .style(BLACK)
            .data( spectrum_iter ),
    )?;


    // ========================= DRAWING THE WEIGHTS =================================
    let mut w = ChartBuilder::on(&weight_root)
        .margin(10)
        .build_cartesian_2d(-50.0..50.0, -50.0..50.0)?;

    // draw triangles
    let weight_iter = (0..triangles.len()).map(|i| w2_sq[(i,i)].sqrt());
    let norm = weight_iter.clone().max_by(|x,y| x.partial_cmp(y).unwrap() ).unwrap() + 0.001;
    for (t, weight) in triangles.iter().zip( weight_iter ) {
        w.draw_series(std::iter::once(Polygon::new(
            t.iter().map(|&i| points[i] ).map(|x| (x[0], x[1]) ).collect::<Vec<_>>(),
            get_color_bwr( weight/norm )
        )))?;
    }

    // draw edges
    let weight_iter = (0..edges.len()).map(|i| w1[(i,i)] );
    let norm = weight_iter.clone().max_by(|x,y| x.partial_cmp(y).unwrap() ).unwrap() + 0.001;
    for (ends, weight) in edges.iter().zip( weight_iter ) {
        w.draw_series(LineSeries::new(
            ends.into_iter().map(|&i| points[i] ).map(|v| (v[0], v[1]) ),
            get_color_bwr( weight/norm ).filled().stroke_width(2)
        ))?;
    }

    // draw vertices
    w.draw_series( 
        points.iter().map(|v| Circle::new((v[0], v[1]), 2, BLACK.filled()))
    )?;



    // ========================= DRAWING THE EIGENVECTORS =================================
    for (i, panel) in eigenvec_panels.iter().enumerate() {
        let mut panel = ChartBuilder::on(&panel)
        .margin(10)
        .build_cartesian_2d(-50.0..50.0, -50.0..50.0)?;

        // draw edges
        let norm = (0..edges.len()).map(|j| eigenvecs[(j,i)].abs() ).max_by(|x,y| x.partial_cmp(y).unwrap() ).unwrap() + 0.001;
        for (j, &edge) in edges.iter().enumerate() {
            panel.draw_series(LineSeries::new(
                edge.into_iter().map(|v| points[v] ).map(|p| (p[0], p[1])),
                get_color_bwr( eigenvecs[(j,i)]/norm ).filled().stroke_width(2)
            ))?;
        }
    
        // draw vertices
        panel.draw_series( 
            points.iter().map(|v| Circle::new((v[0], v[1]), 2, BLACK.filled()))
        )?;
    }
    

    histogram_root.present()?;
    weight_root.present()?;
    eigenvec_root.present()?;

    Ok(())
}

fn get_color_bwr(val: f64) -> RGBColor {
    assert!( -1.0 <= val && val <= 1.0, "val = {val}");
    if val >= 0.0 {
        RGBColor(255, (255.0 * (1.0 - val)) as u8, (255.0 * (1.0 - val)) as u8)
    } else {
        RGBColor((255.0 * (1.0 + val)) as u8, (255.0 * (1.0 + val)) as u8, 255)
    }
}

fn w(e: &[usize; 2], points: &Vec<ConstVector<f64, 2>>)  -> f64 {
    1.0 / (points[e[0]] - points[e[1]]).two_norm()
}