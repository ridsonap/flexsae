# DISCREPANCIES — perbedaan implementasi flexsae vs rumus/dokumentasi publik

Aturan (pola fastsae): **kode adalah sumber kebenaran.** Setiap penyimpangan yang disengaja
dari makalah rujukan atau package pembanding dicatat di sini, di roxygen `@details`, dan di
`NEWS.md`, beserta alasannya. Butir yang belum dipastikan ditandai **TIDAK PASTI**.

Status (28 Okt 2026): `pspline_sae()`, `merf_sae()` dan `benchmark_flexsae()` sudah ada.

## Penyimpangan

| ID | Fungsi | Rujukan | Penyimpangan | Alasan |
|---|---|---|---|---|
| MERF-1 | `merf_sae()` (`R/fh_lm.R`) | Bukhari dkk. (2026), pers. (3)–(5) | σ²ᵥ dari profil likelihood REML (default) atau ML, bukan EM Valdez dkk. | Titik tetap EM = ML (spike S4), jadi `method = "ML"` memberi fit yang sama dengan paper. Default REML agar konsisten dengan `pspline_sae()` dan `sae::eblupFH()`. |
| MERF-2 | `merf_sae()` | Bukhari dkk. (2026), pers. (9), (13) | **Bukan penyimpangan lagi.** Paper eksplisit memakai prediksi OOB; flexsae juga. Area tak-sampel memakai prediksi RF biasa (paper tidak menyebutnya). | Revisi sebelumnya menandai ini TIDAK PASTI; paper sudah dibaca ulang. |
| MERF-3 | `merf_sae(algorithm = "merf")` | Hajjem dkk. (2014); Krennmair & Schmid (2022) | Iterasi RF ↔ FH berofset untuk area-level: adaptasi sendiri, tidak ada di paper mana pun. Kriteria henti hanya pada σ²ᵥ. | Pada `"fhrf"`, σ̂²ᵥ FH linear juga menyerap non-linearitas yang tidak tertangkap (D = 300: 0,430 vs 0,311 setelah iterasi). Eksperimental; dievaluasi di simulasi fase 5. |
| MERF-4 | `merf_sae()` | Bukhari dkk. (2026) | RF memakai variabel RHS formula tanpa transformasi; FH linear memakai suku formula apa adanya. | Paper tidak membahas transformasi; transformasi monoton tidak mengubah RF. |
| PSP-1 | `pspline_sae()` | Ruppert, Wand & Carroll (2003); Wand (2003) | Basis spline berpenalti mgcv (default thin-plate, opsi `"ps"`), bukan basis polinomial terpotong. | Struktur model dan estimand sama; engine REML mgcv (Wood 2011) teruji. |
| BM-1 | `benchmark_flexsae(method = "optimal")` | `saebenchmarking::mse_benchmarking` | Dalam replikat bootstrap, bobot k_d memakai mse_d fit asli (tetap), bukan mse yang dihitung ulang per replikat. | MSE per replikat tidak tersedia tanpa bootstrap bersarang. Estimasi titik sama persis (3e-16). |
| BM-2 | `benchmark_flexsae(method = "optimal")` | `emdi::benchmark(type = "MSE_adj")` | k_d ∝ a_d·mse_d (fastsae, saebenchmarking), bukan k_d ∝ mse_d (emdi). | Mengikuti dua package pembanding lain; perbedaan disebut di Rd. |
| BM-3 | `benchmark_flexsae()` | PLAN §5.4 | Tidak ada MSE linearisasi; hanya bootstrap. | Bootstrap memakai replikat fit-ulang yang sudah ada; untuk difference identik dengan linearisasi bermomen bootstrap. |
| BM-4 | `benchmark_flexsae(national_target = )` | Ghosh & Steorts (2013) | Dua tahap mengikuti algoritma `fastsae::benchmark_sae` (tahap 1 menyesuaikan target grup dengan aturan yang sama), bukan estimator Bayes dua-tahap Ghosh & Steorts. Target nasional dianggap diketahui eksak. | Konsisten dengan fastsae (selisih ≤ 2e-16); paper hanya [V-meta]. |

## Catatan verifikasi

### S1 — FH via mgcv = `sae::eblupFH` (28 Okt 2026)

Data `sae::milk` (43 area), `vardir = SD^2`, REML:

```r
mgcv::gam(yi ~ as.factor(MajorArea) + s(area, bs = "re"),
          weights = 1 / var, scale = 1, method = "REML", data = milk)
```

| Besaran | `sae::eblupFH` | mgcv | Selisih |
|---|---|---|---|
| σ²ᵥ | 0,01855022 | 0,01855033 (`g$scale / g$sp`) | 1,1e-7 |
| β | – | – | maks abs 1,6e-7 |
| EBLUP | – | – | maks relatif 7,8e-7 |

ML (spike yang sama): σ²ᵥ 0,01551752 vs `sae` 0,01551755; EBLUP selisih relatif ≤ 2,6e-7.

Kesimpulan: `weights = 1/ψ` dengan `scale = 1` mengikat varians sampling ke ψ_d, dan
`s(area, bs = "re")` mereproduksi efek area FH. Selisih orde 1e-7 berasal dari kriteria henti
optimizer. Akan dijadikan test regresi (toleransi 1e-5 untuk σ²ᵥ, 1e-6 relatif untuk EBLUP).

### S2 — `emdi::fh` REML berbeda ~1,5e-4 relatif (28 Okt 2026)

Pada `milk`, `emdi::fh(method = "reml")` memberi σ²ᵥ = 0,01854753, sedangkan `sae::eblupFH`,
`fastsae::eblup_fh` dan `pspline_sae()` memberi 0,0185502–0,0185503. Penyebab (dibaca dari kode
emdi 2.2.3, fungsi internal `Reml`): emdi memaksimalkan likelihood REML dengan `stats::optimize()`
pada toleransi default (`.Machine$double.eps^0.25` ≈ 1,2e-4); argumen `tol` emdi tidak berlaku untuk
FH REML tanpa korelasi spasial (hasil sama untuk `tol = 1e-4` dan `1e-8`). Jadi selisihnya adalah
presisi optimizer emdi, bukan perbedaan model. Test terhadap emdi memakai toleransi 5e-4 (σ²ᵥ) dan
1e-4 (estimasi area); terhadap sae dan fastsae tetap 1e-5 / 1e-6.

### S3 — refit bootstrap tidak memakai `mgcv::gam(G = )` (28 Okt 2026)

Mengganti `G$y` pada objek `gam(..., fit = FALSE)` lalu memanggil `gam(G = G)` memberi parameter
penghalus berbeda dari fit penuh pada respons yang sama (sp 46,7 vs 56,5 pada `milk`). Bootstrap
`pspline_sae()` karena itu mem-fit ulang penuh dari data. Penyebab di mgcv belum ditelusuri.

### S4 — titik tetap EM Bukhari dkk. (2026) = ML (28 Okt 2026)

EM pers. (3)–(5) pada `sae::milk`, dengan langkah β GLS pada y − v̂ (pers. (4) di paper tercetak sama
dengan (3), jadi langkah β tidak tertulis): 77 iterasi, σ²ᵥ = 0,01551751 vs `sae::eblupFH(ML)`
0,01551755; |Δβ| maks 7,5e-8. Kesimpulan: FH–EM paper adalah estimator ML. Skrip spike ada di
scratchpad sesi; hasilnya dijadikan dasar MERF-1.

### S5 — pemecah FH satu komponen varians vs mgcv (28 Okt 2026)

`.fit_fh_lm()` (grid 101 titik di [0, 5·var(y)], lalu `optimize()` dengan tol 1e-10·batas atas) pada
`milk`: REML σ²ᵥ 0,01855033 (sae 0,01855022), ML 0,01551751 (sae 0,01551755); |Δβ| ≤ 1,6e-7, |ΔSE| ≤ 1,7e-7.
Waktu per fit: 0,007 dtk (milk); 0,006 dtk pada D = 300 (rata-rata 20 fit) vs mgcv 2,57 dtk, σ²ᵥ identik (0,4295193). Tanpa efek tetap,
REML = ML dan sama dengan maksimisasi langsung (test). Tabel koefisien memakai uji z (konvensi `sae`).

### S6 — pilot kalibrasi bootstrap MSE `merf_sae(algorithm = "fhrf")` (28 Okt 2026)

Skrip `benchmarks/pilot/merf_mse_pilot.R`: x tetap, D = 100, θ_d = 2 + sin(1,5 x1) + 0,5 x2 + v_d
(σᵥ = 0,5), ψ_d ~ U(0,2; 0,6), R = 100 replikasi, B = 50, `num.trees = 200`, REML.

| Varian | MSE empiris (rata-rata area) | Bootstrap (rata-rata) | Rasio | RB area: median [P10, P90] |
|---|---|---|---|---|
| replace | 0,1988 | 0,2299 | 1,156 | 0,159 [−0,022; 0,472] |
| correct | 0,1957 | 0,2311 | 1,181 | 0,170 [−0,017; 0,424] |

Bootstrap **melebihkan** MSE ~16–18% (konservatif). **Dugaan penyebab (belum diuji):** σ̂²ᵥ dari FH
linear menyerap non-linearitas (±0,43 pada data sejenis vs nilai benar 0,25), sehingga θ* = m̂ + v*
terlalu menyebar. Ini sifat skema Bukhari dkk. pers. (15)–(17) bersama σ̂²ᵥ linear, bukan bug yang
teridentifikasi. Yang perlu diuji di fase 5 (ADEMP dengan MCSE): (a) `algorithm = "merf"` (σ²ᵥ
diestimasi ulang), (b) DGP linear (rasio seharusnya ≈ 1), (c) D lebih besar. Pilot ini tidak
melaporkan MCSE dan tidak boleh dikutip sebagai hasil simulasi.

### S7 — `benchmark_flexsae()` vs fastsae, emdi, saebenchmarking (28 Okt 2026)

Aturan (vektor estimasi, 30 area, 3 grup; skrip di test `test_benchmark.R`):

| Pembanding | Kasus | Selisih maks |
|---|---|---|
| `fastsae::benchmark_sae` | ratio, difference × mean/total × satu/banyak grup × dua tahap (dengan/tanpa target grup) | ≤ 2,2e-16 |
| idem | logit | ≤ 2,2e-11 (fastsae `uniroot` tol 1e-9, flexsae 1e-12) |
| `fastsae:::.benchmark_worker` | optimal dengan mse | 0 |
| `saebenchmarking::sae_benchmarking` | optimum (phi_source = "mse") | 3,3e-16 |
| `emdi::benchmark` | raking = difference, ratio | 2,2e-16 |

MSE (milk, bobot N ~ Poisson(300), `pspline_sae` linear, B = 2000, target `"direct"`):

| Metode | Rata-rata kenaikan MSE (ours, bootstrap) | Pembanding | MSE terbenchmark ours / saebenchmarking per area |
|---|---|---|---|
| difference | 0,000078 | g₄ Steorts & Ghosh (2013) via `saebenchmarking`: 0,000077 | median 0,972 [0,908; 1,027] |
| ratio | 0,000084 | bootstrap `saebenchmarking` (B = 1000) | median 0,989 [0,912; 1,122] |

Selisih level (~3%) berasal dari MSE dasar bootstrap vs Prasad-Rao (rasio 0,98 di fase 1), bukan dari benchmark.
