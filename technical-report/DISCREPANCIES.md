# DISCREPANCIES — perbedaan implementasi flexsae vs rumus/dokumentasi publik

Aturan (pola fastsae): **kode adalah sumber kebenaran.** Setiap penyimpangan yang disengaja
dari makalah rujukan atau package pembanding dicatat di sini, di roxygen `@details`, dan di
`NEWS.md`, beserta alasannya. Butir yang belum dipastikan ditandai **TIDAK PASTI**.

Status: belum ada kode metode (28 Okt 2026). Butir di bawah adalah penyimpangan yang sudah
**direncanakan** di `PLAN.md`; nomor baris ditambahkan saat kodenya ditulis.

## Direncanakan

| ID | Fungsi | Rujukan | Penyimpangan | Alasan |
|---|---|---|---|---|
| MERF-1 | `merf_sae()` | Bukhari dkk. (2026), pers. (3)–(5) | σ²ᵥ diestimasi dengan REML Fay-Herriot dengan offset f̂, bukan EM Valdez dkk. | Rujukan EM belum diverifikasi; REML FH via mgcv sudah terbukti sama dengan `sae::eblupFH` (catatan S1). |
| MERF-2 | `merf_sae()` | Bukhari dkk. (2026), pers. (11)–(14) | Langkah efek acak memakai prediksi RF out-of-bag, bukan in-sample. **TIDAK PASTI** apakah Bukhari dkk. memakai OOB. | Prediksi in-sample RF menyerap galat sampling y_d sehingga σ²ᵥ terlalu kecil; OOB mengikuti Krennmair & Schmid (2022). |
| PSP-1 | `pspline_sae()` | Ruppert, Wand & Carroll (2003); Wand (2003) | Basis spline berpenalti mgcv (default thin-plate, opsi `"ps"`), bukan basis polinomial terpotong. | Struktur model dan estimand sama; engine REML mgcv (Wood 2011) teruji. |

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
