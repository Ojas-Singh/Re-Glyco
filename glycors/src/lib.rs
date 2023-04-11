use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use std::f64::consts::PI;
use rand::distributions::{Uniform, Distribution};
use rayon::prelude::*;
use rand::Rng;

#[inline(always)]
fn length(v: &Vec<f64>) -> f64 {
    let sum = v.iter().map(|x| x * x).sum::<f64>();
    sum.sqrt()
}

#[inline(always)]
fn matrix_multiply(a: &Vec<Vec<f64>>, b: &Vec<f64>) -> Vec<f64> {
    assert_eq!(b.len(), 3);

    let mut result = vec![0.0; 3];

    for i in 0..3 {
        for j in 0..3 {
            result[i] += a[i][j] * b[j];
        }
    }

    result
}

#[inline(always)]
fn normalize(v: &Vec<f64>) -> Vec<f64> {
    let len = length(v);
    vec![v[0] / len, v[1] / len, v[2] / len]
}

#[inline(always)]
fn dot(v: &Vec<f64>, w: &Vec<f64>) -> f64 {
    v.iter().zip(w.iter()).map(|(x, y)| x * y).sum()
}

#[inline(always)]
fn cross(v: &Vec<f64>, w: &Vec<f64>) -> Vec<f64> {
    vec![
        v[1] * w[2] - v[2] * w[1],
        v[2] * w[0] - v[0] * w[2],
        v[0] * w[1] - v[1] * w[0],
    ]
}

#[inline(always)]
fn subtract(a: &Vec<f64>, b: &Vec<f64>) -> Vec<f64> {
    vec![a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline(always)]
fn norm(v: &Vec<f64>) -> f64 {
    let sum = v.iter().map(|x| x * x).sum::<f64>();
    sum.sqrt()
}

#[inline(always)]
fn cdist(a: &Vec<Vec<f64>>, b: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    assert_eq!(a[0].len(), b[0].len());

    let (a_rows, a_cols) = (a.len(), a[0].len());
    let (b_rows, _) = (b.len(), b[0].len());

    let c: Vec<Vec<f64>> = (0..a_rows)
        .into_par_iter()
        .map(|i| {
            let mut c_row = vec![0.0; b_rows];
            for j in 0..b_rows {
                let mut acc = 0.0;
                for k in 0..a_cols {
                    let diff = a[i][k] - b[j][k];
                    acc += diff * diff;
                }
                c_row[j] = acc.sqrt();
            }
            c_row
        })
        .collect();

    c
}

#[inline(always)]
fn steric_fast(garr: &Vec<Vec<f64>>, parr: &Vec<Vec<f64>>) -> f64 {
    let c = cdist(garr, parr);
    let threshold = 1.7;

    let mut result = 1.0;

    for i in 3..c.len() {
        for elem in &c[i] {
            if *elem < threshold {
                result += 200.0 * (-elem.powi(2)).exp();
            }
        }
    }

    result
}

#[inline(always)]
fn rr(
    phi: f64,
    psi: f64,
    cb: usize,
    cg: usize,
    nd2: usize,
    c1: usize,
    o5: usize,
    garr: &mut Vec<Vec<f64>>,
    parr: &Vec<Vec<f64>>,
) {
    let rad_phi = PI / 180.0 * (phi - fastest_dihedral(&parr[cb], &parr[cg], &parr[nd2], &garr[c1]));
    let m1 = rotation_matrix(&subtract(&parr[nd2], &parr[cg]), rad_phi);

    for row in garr.iter_mut() {
        for i in 0..row.len() {
            row[i] -= parr[nd2][i];
        }
    }
    for row in garr.iter_mut() {
        let new_row = matrix_multiply(&m1,&row);
        *row = new_row;
    }
    for row in garr.iter_mut() {
        for i in 0..row.len() {
            row[i] += parr[nd2][i];
        }
    }

    let rad_psi = PI / 180.0 * (psi - fastest_dihedral(&parr[cg], &parr[nd2], &garr[c1], &garr[o5]));
    let m2 = rotation_matrix(&subtract(&garr[c1], &parr[nd2]), rad_psi);

    for row in garr.iter_mut() {
        for i in 0..row.len() {
            row[i] -= parr[nd2][i];
        }
    }
    for row in garr.iter_mut() {
        let new_row = matrix_multiply(&m2,&row);
        *row = new_row;
    }
    
    for row in garr.iter_mut() {
        for i in 0..row.len() {
            row[i] += parr[nd2][i];
        }
    }
}

#[pyfunction]
pub fn opt_genetic(
    cb: usize,
    cg: usize,
    nd2: usize,
    c1: usize,
    o5: usize,
    garr: Vec<Vec<f64>>,
    parr: Vec<Vec<f64>>,
    phisd: (f64, f64),
    psisd: (f64, f64),
) -> PyResult<(f64, f64, f64)> {
    let population_size = 100;
    let generations = 20;
    let mutation_rate = 0.1;
    let mut rng = rand::thread_rng();
    if phisd.0 >= phisd.1 || psisd.0 >= psisd.1 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid input ranges for phisd and/or psisd"));
    }
    

    let phi_range = Uniform::new(phisd.0, phisd.1);
    let psi_range = Uniform::new(psisd.0, psisd.1);

    // Generate initial population
    let mut population: Vec<(f64, f64)> = (0..population_size)
        .map(|_| (phi_range.sample(&mut rng), psi_range.sample(&mut rng)))
        .collect();

    let mut best_individual = (0.0, 0.0);
    let mut best_fitness = 1e15;

    for _ in 0..generations {
        // Evaluate fitness of the population
        let fitness: Vec<f64> = population
            .iter()
            .map(|&(phi, psi)| {
                let mut garr_temp = garr.clone();
                rr(phi, psi, cb, cg, nd2, c1, o5, &mut garr_temp, &parr);
                steric_fast(&garr_temp, &parr)
            })
            .collect();

        // Find the best individual in the current population
        if let Some((idx, &fit)) = fitness.iter().enumerate().min_by(|a, b| a.1.partial_cmp(b.1).unwrap()) {
            if fit < best_fitness {
                best_fitness = fit;
                best_individual = population[idx];
            }
        }

        // Select parents and create offspring
        let mut offspring: Vec<(f64, f64)> = vec![];

        while offspring.len() < population_size {
            let parents = select_parents(&population, &fitness, &mut rng);
            let children = crossover(parents.0, parents.1, &mut rng);

            offspring.push(mutate(children.0, mutation_rate, &phi_range, &psi_range, &mut rng));
            offspring.push(mutate(children.1, mutation_rate, &phi_range, &psi_range, &mut rng));
        }

        population = offspring;
    }

    Ok((best_individual.0, best_individual.1, best_fitness))
}


fn select_parents<'a>(
    population: &'a [(f64, f64)],
    fitness: &[f64],
    rng: &mut impl Rng,
) -> (&'a (f64, f64), &'a (f64, f64)) {
    let total_fitness: f64 = fitness.iter().sum();
    let mut selection_probabilities: Vec<f64> = fitness.iter().map(|&fit| total_fitness - fit).collect();

    for i in 1..selection_probabilities.len() {
        selection_probabilities[i] += selection_probabilities[i - 1];
    }

    let parent_idx1 = selection_probabilities
        .iter()
        .position(|&x| x >= rng.gen_range(0.0..*selection_probabilities.last().unwrap()))
        .unwrap();
        let parent_idx2 = selection_probabilities
        .iter()
        .position(|&x| x >= rng.gen_range(0.0..*selection_probabilities.last().unwrap()))
        .unwrap();
        (&population[parent_idx1], &population[parent_idx2])
    }

fn crossover(
    parent1: &(f64, f64),
    parent2: &(f64, f64),
    rng: &mut impl Rng,
    ) -> ((f64, f64), (f64, f64)) {
    let crossover_point = rng.gen_range(0..=1);
    if crossover_point == 0 {
        ((parent1.0, parent2.1), (parent2.0, parent1.1))
    } else {
        ((parent1.0, parent1.1), (parent2.0, parent2.1))
    }
    }

fn mutate(
        individual: (f64, f64),
        mutation_rate: f64,
        phi_range: &Uniform<f64>,
        psi_range: &Uniform<f64>,
        rng: &mut impl Rng,
        ) -> (f64, f64) {
        let mut mutated_individual = individual;
        if rng.gen::<f64>() < mutation_rate {
            mutated_individual.0 = phi_range.sample(rng);
        }
        
        if rng.gen::<f64>() < mutation_rate {
            mutated_individual.1 = psi_range.sample(rng);
        }
        
        mutated_individual
    }


#[pyfunction]
pub fn opt(
    cb: usize,
    cg: usize,
    nd2: usize,
    c1: usize,
    o5: usize,
    garr: Vec<Vec<f64>>,
    parr: Vec<Vec<f64>>,
    phisd: (f64, f64),
    psisd: (f64, f64),
) -> PyResult<(f64, f64, f64)> {
    let mut phif = 0.0;
    let mut psif = 0.0;
    let mut r = 1e15;
    let mut rng = rand::thread_rng();
    let phi_range = Uniform::new(phisd.0, phisd.1);
    let psi_range = Uniform::new(psisd.0, psisd.1);

    for _ in 0..1000 {
        let phi = phi_range.sample(&mut rng);
        let psi = psi_range.sample(&mut rng);
        let mut garr_temp = garr.clone();
        rr(phi, psi, cb, cg, nd2, c1, o5, &mut garr_temp, &parr);
        let ri = steric_fast(&garr_temp, &parr);
        if ri < r {
            phif = phi;
            psif = psi;
            r = ri;
        }
    }

    Ok((phif, psif, r))
}


fn rotation_matrix(axis: &Vec<f64>, theta: f64) -> Vec<Vec<f64>> {
    let axis_norm = {
        let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        [axis[0] / len, axis[1] / len, axis[2] / len]
    };
    let (a, b, c, d) = (
        (theta / 2.0).cos(),
        -axis_norm[0] * (theta / 2.0).sin(),
        -axis_norm[1] * (theta / 2.0).sin(),
        -axis_norm[2] * (theta / 2.0).sin(),
    );
    let (aa, bb, cc, dd) = (a * a, b * b, c * c, d * d);
    let (bc, ad, ac, ab, bd, cd) = (b * c, a * d, a * c, a * b, b * d, c * d);
    vec![
        vec![aa + bb - cc - dd, 2.0 * (bc + ad), 2.0 * (bd - ac)],
        vec![2.0 * (bc - ad), aa + cc - bb - dd, 2.0 * (cd + ab)],
        vec![2.0 * (bd + ac), 2.0 * (cd - ab), aa + dd - bb - cc],
    ]
}


fn fastest_dihedral(p0: &Vec<f64>, p1: &Vec<f64>, p2: &Vec<f64>, p3: &Vec<f64>) -> f64 {
    let b1 = vec![
        p2[0] - p1[0],
        p2[1] - p1[1],
        p2[2] - p1[2],
    ];
    let b1_norm = normalize(&b1);
    let b0 = vec![
        -(p1[0] - p0[0]),
        -(p1[1] - p0[1]),
        -(p1[2] - p0[2]),
    ];
    let b2 = vec![
        p3[0] - p2[0],
        p3[1] - p2[1],
        p3[2] - p2[2],
    ];
    let v = vec![
        b0[0] - dot(&b0, &b1_norm) * b1_norm[0],
        b0[1] - dot(&b0, &b1_norm) * b1_norm[1],
        b0[2] - dot(&b0, &b1_norm) * b1_norm[2],
    ];
    let w = vec![
        b2[0] - dot(&b2, &b1_norm) * b1_norm[0],
        b2[1] - dot(&b2, &b1_norm) * b1_norm[1],
        b2[2] - dot(&b2, &b1_norm) * b1_norm[2],
    ];
    let x = dot(&v, &w);
    let y = dot(&cross(&b1_norm, &v), &w);
    180.0 * y.atan2(x) / PI
}


fn fastest_angle(p0: &Vec<f64>, p1: &Vec<f64>, p2: &Vec<f64>) -> f64 {
    let v0 = subtract(&p0, &p1);
    let v1 = subtract(&p2, &p1);
    let cosine_angle = dot(&v0, &v1) / (norm(&v0) * norm(&v1));
    let angle = cosine_angle.acos();
    angle.to_degrees()
}


#[pymodule]
fn glycors(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(opt))?; // Fix: Remove the '?' inside the parentheses
    m.add_wrapped(wrap_pyfunction!(opt_genetic))?;
    Ok(())
}


