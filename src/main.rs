use std::mem;
extern crate crossbeam;
extern crate scoped_pool;

fn main() {
    let residues = [1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
                    71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 121, 127, 131,
                    137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181, 187,
                    191, 193, 197, 199, 209, 211];
    let val = 2_000_000_000;
    let md = 210;
    let rescnt = 48;

    println!("val = {}, mod = {}, rescnt = {}", val, md, rescnt);

    let mut posn = [0; 210];
    for i in 1..rescnt {
        posn[residues[i]] = i-1;
    }
    posn[1] = rescnt-1;

    let mut modk; let mut r; let mut k;

    let num = val-1 | 1;           // if value even number then subtract 1
    k = num/md; modk = md*k; r=1;  // kth residue group, base num value
    // find last pc position <= num
    while num >= modk+residues[r] {
        r += 1;
    }
    let maxpcs  = k*rescnt + r-1; // maximum number of prime candidates

    let prms: Vec<u8> = vec![0; maxpcs];

    println!("num = {}, k = {}, modk = {}, maxpcs = {}", num, k, modk, maxpcs);

    let sqrt_n = (num as f32).sqrt() as usize;

    modk=0; r=0; k=0;

    let pool = scoped_pool::Pool::new(3);
    // sieve to eliminate nonprimes from primes prms array
    for i in 0..maxpcs {
        r += 1; if r > rescnt {r = 1; modk += md; k += 1;};
        if prms[i] == 1 {
            continue;
        }
        let prm_r = residues[r];
        let prime = modk + prm_r;
        if prime > sqrt_n {
            break;
        }
        let prmstep = prime * rescnt;
        pool.scoped(|scope| {
            for ri in &residues[1..rescnt+1] {
                let prms = &prms;
                scope.execute(move || {
                    let raw_prms = unsafe { mem::transmute::<*const u8, *mut u8>(prms.as_ptr()) };
                    let prod = prm_r * ri;
                    let mut j = (k*(prime + ri) + (prod-2)/md)*rescnt + posn[prod % md];
                    while j < maxpcs {
                        unsafe { *raw_prms.offset(j as isize) = 1; }
                        j += prmstep;
                    }
                });
            }
        });
    }
    // the prms array now has all the positions for primes r1..N
    // extract prime numbers and count from prms into prims array
    let mut prmcnt = 4;
    modk=0; r=0;
    for i in 0..maxpcs {
        r += 1;
        if r > rescnt {
            r = 1;
            modk += md;
        }
        if prms[i] == 0 {
            prmcnt += 1;
        }
    }
    println!("{}", prmcnt);
}
