use rand::Rng;
use std::collections::HashSet;

fn power_mod(base_val: u64, exponent: u64, modulus: u64) -> u64 {
    let mut result_val = 1;
    let mut current_base = base_val % modulus;
    let mut current_exp = exponent;

    while current_exp > 0 {
        if current_exp & 1 == 1 {
            result_val = (result_val * current_base) % modulus;
        }
        current_base = (current_base * current_base) % modulus;
        current_exp >>= 1;
    }
    result_val
}

fn compute_poly_value(coeffs: &[u64], x_val: u64, field_mod: u64) -> u64 {
    coeffs.iter().rev().fold(0, |accumulator, &coeff| {
        (accumulator * x_val + coeff) % field_mod
    })
}

fn create_poly_commitments(poly_coeffs: &[u64], generator: u64, prime_mod: u64) -> Vec<u64> {
    poly_coeffs
        .iter()
        .map(|&coeff| power_mod(generator, coeff, prime_mod))
        .collect()
}

fn validate_share(
    x_val: u64,
    y_val: u64,
    commitments: &[u64],
    generator: u64,
    prime_mod: u64,
) -> bool {
    let left_side = power_mod(generator, y_val, prime_mod);

    let right_side = commitments
        .iter()
        .enumerate()
        .fold(1, |accumulator, (idx, &commit_val)| {
            let exponent = x_val.pow(idx as u32);
            (accumulator * power_mod(commit_val, exponent, prime_mod)) % prime_mod
        });

    left_side == right_side
}

fn recover_secret_value(share_points: &[(u64, u64)], field_mod: u64) -> u64 {
    let mut reconstructed = 0;

    for (i, &(x_i, y_i)) in share_points.iter().enumerate() {
        let mut lagrange_coeff = 1;

        for (j, &(x_j, _)) in share_points.iter().enumerate() {
            if i != j {
                let numerator = field_mod + x_j;
                let denominator = field_mod + x_j - x_i;
                let inv_denominator = modular_inverse(denominator % field_mod, field_mod);

                lagrange_coeff =
                    lagrange_coeff * numerator % field_mod * inv_denominator % field_mod;
            }
        }
        reconstructed = (reconstructed + y_i * lagrange_coeff % field_mod) % field_mod;
    }

    reconstructed
}

fn modular_inverse(value: u64, modulus: u64) -> u64 {
    let mut gcd_pair = (modulus, value);
    let mut coeff_pair = (0, 1);

    while gcd_pair.1 != 0 {
        let quotient = gcd_pair.0 / gcd_pair.1;
        gcd_pair = (gcd_pair.1, gcd_pair.0 % gcd_pair.1);
        coeff_pair = (
            coeff_pair.1,
            coeff_pair.0 - (quotient as i64) * coeff_pair.1,
        );
    }

    ((coeff_pair.0 + modulus as i64) % modulus as i64) as u64
}

fn main() {
    let prime_modulus = 23;
    let cyclic_generator = 2;
    let field_size = 11;
    let original_secret = 2;

    println!("Original secret value: {}", original_secret);

    let total_participants = 7;
    let minimum_shares = 3;
    println!("Scheme: {}-out-of-{}", minimum_shares, total_participants);

    let mut random_gen = rand::rng();
    let mut poly_coeffs: Vec<u64> = (0..minimum_shares - 1)
        .map(|_| random_gen.random_range(1..field_size))
        .collect();

    poly_coeffs.insert(0, original_secret);

    println!(
        "Generated polynomial: {} + {}x + {}x²",
        poly_coeffs[0], poly_coeffs[1], poly_coeffs[2]
    );

    let poly_commitments = create_poly_commitments(&poly_coeffs, cyclic_generator, prime_modulus);
    println!("Polynomial commitments: {:?}", poly_commitments);

    let distributed_shares: Vec<_> = (1..=total_participants)
        .map(|participant_id| {
            (
                participant_id as u64,
                compute_poly_value(&poly_coeffs, participant_id as u64, field_size),
            )
        })
        .collect();

    println!("Distributed shares: {:?}", distributed_shares);

    let mut chosen_indices = HashSet::new();
    while chosen_indices.len() < minimum_shares {
        chosen_indices.insert(random_gen.random_range(0..total_participants));
    }

    let selected_shares: Vec<_> = chosen_indices
        .iter()
        .map(|&idx| distributed_shares[idx])
        .collect();

    println!("Selected shares for reconstruction: {:?}", selected_shares);

    for &(x_coord, y_coord) in &selected_shares {
        if validate_share(
            x_coord,
            y_coord,
            &poly_commitments,
            cyclic_generator,
            prime_modulus,
        ) {
            println!("✓ Share ({}, {}) verified successfully", x_coord, y_coord);
        } else {
            println!("✗ Share ({}, {}) failed verification!", x_coord, y_coord);
        }
    }

    let recovered_secret = recover_secret_value(&selected_shares, field_size);
    println!("Recovered secret value: {}", recovered_secret);

    assert_eq!(original_secret, recovered_secret);
    println!("✓ Secret reconstruction successful!");
}
