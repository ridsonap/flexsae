# PLAN — package `flexsae`: Small Area Estimation area-level dengan hubungan non-linear

Status: **rencana, belum ada kode metode**. Revisi 28 Okt 2026: lingkup v0.1 **area-level saja**.
Rujukan: `technical-report/references.bib` (68 entri, semuanya terverifikasi) dan
`technical-report/TODO_REFERENCES.md` (yang belum terverifikasi, tidak boleh dikutip).
Penyimpangan yang disengaja: `technical-report/DISCREPANCIES.md`.

Penanda: **[V]** metadata terverifikasi dan isi yang dipakai sudah dibaca; **[V-meta]** metadata
terverifikasi, isi detail belum dibaca (jangan implementasikan persamaannya sebelum dibaca);
**[BELUM]** tidak terverifikasi, tidak dikutip; **[dugaan]** penalaran sendiri yang harus
dibuktikan dengan tes atau derivasi.

> Catatan: prompt awal terpotong di "Sensitivitas: ... model salah spesifikasi (linear".
> Butir sensitivitas di §5.6 sebagian adalah tebakan; perbarui bila sisa prompt tersedia.

---

## 0. Keputusan

| Butir | Keputusan |
|---|---|
| Nama | `flexsae`, tersedia di CRAN (dicek 28 Okt 2026). Repo publik `github.com/ridsonap/flexsae`. |
| Lingkup v0.1 | **Area-level**: data per area berupa estimasi langsung $y_d$, varians sampling $\psi_d$, kovariat $x_d$. |
| Metode | **`pspline_sae()`** (Fay-Herriot semiparametrik P-spline) dan **`merf_sae()`** (Fay-Herriot + random forest). BART-FH hanya pembanding lewat `nlfh`. |
| Penamaan | Fungsi model `*_sae`; utilitas `*_flexsae` (`benchmark_flexsae()`, `select_flexsae()`, `compare_flexsae()`, `importance_flexsae()`, `partial_flexsae()`, `sim_flexsae_data()`) agar tidak menutupi `benchmark_sae()`, `step_sae()`, `compare_sae()`, `map_sae()` milik fastsae. |
| Engine | `mgcv` (pspline, dan langkah FH di merf), `ranger` (merf). R murni dulu. |
| MSE | Bootstrap parametrik yang mengulang seluruh fit untuk kedua metode. |
| Benchmark | `benchmark_flexsae()` post-hoc + MSE terbenchmark (linearisasi dan bootstrap ulang penuh); self-benchmark lewat augmentasi untuk kedua metode (merf: [dugaan]). |
| Pembanding linear | `sae` 1.3 (`eblupFH`, `mseFH`), `emdi` 2.2.3 (`fh`, `benchmark`), `fastsae` (`eblup_fh`, `benchmark_sae`, `step_sae`), `saebenchmarking` (MSE terbenchmark FH), `nlfh` (BART-FH). |

---

## 1. Acuan dari fastsae (`origin/main` 6bdc2a6)

Ditiru: satu kelas hasil dengan tabel berkolom baku yang dibaca semua metode S3 (`.sae_table()`
di `R/utils.R`); pemisahan tegas self-benchmark (argumen fungsi model) vs post-hoc
(`R/benchmark.R`, kelas `fastsae_benchmark`); MSE terbenchmark dengan operator matriks MSE penuh
(`.mse_operator()`); `step_sae()` (IC dari log-likelihood penuh, fit ML selama pencarian, refit
dengan metode asal, validasi simulasi di `benchmarks/validation/step_sae.R`); `autoplot()` dengan
`type =` dan helper `.autoplot_*`; pengukuran peak RSS satu fit per proses
(`benchmarks/eval/bench_resources.R`, `ps`, minus baseline); `technical-report/` dengan bib
terverifikasi, `TODO_REFERENCES.md`, `DISCREPANCIES.md`.

Tidak diwarisi begitu saja: MSE terbenchmark fastsae menganggap target pasti (tipe (i) tidak
tertangani), dan self-benchmark fastsae menormalkan bobot atas domain tersampel sementara post-hoc
atas semua domain. Jebakan yang dihindari: `Collate:` dan `NAMESPACE` tulis tangan (flexsae memakai
roxygen untuk `NAMESPACE`, penyimpangan sengaja); skrip benchmark harus memakai paket terinstal.

---

## 2. Metode

### 2.0 Literatur area-level non-linear

| Rujukan | Isi | Status |
|---|---|---|
| Fay & Herriot (1979) | Model area-level dasar | [V-meta] |
| Ruppert, Wand & Carroll (2003); Wand (2003) | Representasi mixed model untuk spline berpenalti | [V-meta] |
| Wood (2011, 2017) | REML untuk GAM, mgcv | [V-meta] |
| Giusti, Marchetti, Pratesi & Salvati (2012) | FH P-spline frekuentis | **[BELUM]** (hanya sitasi sekunder) |
| **Bukhari, Notodiputro, Indahwati & Fitrianto (2026)** *Ing. Syst. d'Inf.* 31(6):1811–1823 | FH + random forest (FHRF-C, FHRF-R), bootstrap MSE | [V] (halaman penerbit, persamaan (1)–(17)) |
| Bukhari dkk. (2025) *BIO Web Conf.* 186:02012 | FH termodifikasi dengan ML (prosiding) | [V-meta] |
| Hajjem, Bellavance & Larocque (2014); Krennmair & Schmid (2022) | MERF (unit-level), dasar algoritma | [V-meta] / [V] (arXiv v2) |
| Parker (2024) JOS 40(2) | FH non-linear dengan random-weight NN (Bayes) | [V] (persamaan model) |
| Parker & Eideh (2026) JSSAM | BART-FH | [V-meta] |
| Sugasawa, Kubokawa & Rao (2019); Polettini (2017) | HB dengan link tak diketahui; FH semiparametrik Bayes | [V-meta] |

Catatan jujur: rujukan frekuentis FH P-spline belum terverifikasi, jadi `pspline_sae()` dirujuk
ke dua komponen yang terverifikasi (FH + representasi mixed model P-spline). Rujukan
`merf_sae()` (Bukhari dkk. 2026) terbit di venue berprofil rendah untuk metodologi statistik;
dasar algoritmanya tetap MERF Hajjem dkk. (2014) dan Krennmair & Schmid (2022), dan Rd harus
menyebut keduanya.

### 2.1 Notasi dan data

Area $d = 1,\dots,D$; estimasi langsung $y_d$, varians sampling $\psi_d > 0$ dianggap diketahui
(boleh hasil GVF), kovariat area $x_d$, ukuran populasi $N_d$ (untuk benchmark). Area tak-sampel:
$y_d$ = NA, $x_d$ tersedia. Estimand: $\theta_d$.

### 2.2 `pspline_sae()` — Fay-Herriot semiparametrik P-spline

- **Model.** $y_d = \theta_d + e_d$, $e_d \sim N(0,\psi_d)$; $\theta_d = \beta_0 + \sum_k f_k(x_{dk}) + x_d^{L\prime}\beta + v_d$,
  $v_d \sim N(0,\sigma_v^2)$; $f_k$ spline berpenalti, penalti ekuivalen dengan efek acak (Ruppert dkk. 2003;
  Wand 2003). Kasus tanpa smooth = FH linear.
- **Fit.** `mgcv::gam(y ~ s(x1) + x2 + s(domain, bs = "re"), weights = 1/psi, scale = 1, method = "REML")`.
  **Spike S1 sudah lulus**: pada `sae::milk`, hasilnya sama dengan `sae::eblupFH(REML)` (σ²ᵥ selisih
  1,1e-7, EBLUP selisih relatif ≤ 7,8e-7; σ²ᵥ = `g$scale / g$sp`). Detail: `DISCREPANCIES.md` §S1.
  `method = "ML"` untuk seleksi.
- **Prediktor.** $\hat\theta_d = \hat\beta_0 + \sum_k\hat f_k(x_{dk}) + x_d^{L\prime}\hat\beta + \hat v_d$; tak-sampel $\hat v_d = 0$
  (prediksi dengan `predict(..., exclude = "s(domain)")`).
- **MSE.** Bootstrap parametrik naif (Hall & Maiti 2006 [V-meta, abstrak]; skema yang sama dengan
  Bukhari dkk. 2026 pers. 15–17 [V]): $\theta^*_d = \hat m(x_d) + v^*_d$, $v^*_d \sim N(0,\hat\sigma_v^2)$, $y^*_d = \theta^*_d + e^*_d$,
  fit ulang **termasuk parameter penghalus**, $\widehat{\mathrm{MSE}}_d = B^{-1}\sum_b(\hat\theta^*_d - \theta^*_d)^2$. Koreksi bias
  (double bootstrap) Hall & Maiti belum diimplementasikan. Analitik multi-komponen ditunda.
  *Koreksi 28 Okt 2026:* revisi sebelumnya menyebut González-Manteiga dkk. (2008, JSCS) — makalah itu
  untuk model nested-error unit-level, bukan FH; Butar & Lahiri (2003) mengusulkan ukuran bootstrap
  terkoreksi bias untuk model campuran umum, bukan skema naif ini. Keduanya tidak dikutip untuk `pspline_sae()`.
  ✅ Fase 1: untuk FH linear, rata-rata MSE bootstrap / Prasad–Rao (`sae::mseFH`) = 0,98 pada `milk`, B = 400
  (per area 0,85–1,16).
- **Asumsi/risiko.** Normalitas $v_d$ dan $e_d$; $\psi_d$ benar; satu observasi per area membatasi
  fleksibilitas: basis `k` dibatasi menurut $D$ dan diuji di simulasi (perkiraan: kalah dari FH
  linear bila $D$ kecil dan $f$ linear).
- **Penyimpangan.** Basis mgcv, bukan polinomial terpotong (DISCREPANCIES PSP-1).

### 2.3 `merf_sae()` — Fay-Herriot dengan random forest

Sumber: Bukhari dkk. (2026) [V], Bagian 2.3–2.6, pers. (6)–(17) (dibaca ulang 28 Okt 2026); algoritma
MERF Hajjem dkk. (2014), Krennmair & Schmid (2022).

> **Koreksi atas revisi PLAN sebelumnya.** Revisi sebelumnya menulis algoritma Bukhari dkk. sebagai
> iterasi RF ↔ FH sampai konvergen. Itu keliru. Paper menjalankan **satu putaran**: FH–EM linear →
> RF sekali (prediksi OOB) → EBLUP dengan σ̂²ᵥ dari FH linear. Paper juga eksplisit memakai OOB
> ("rely solely on OOB estimates", pers. 9), jadi butir TIDAK PASTI di MERF-2 terjawab.

- **Varian** `variant = c("replace", "correct")`:
  - `"replace"` (FHRF-R, pers. 11–14): RF pada $y_d - \hat v_d$; $\hat m_d = \hat f(x_d)$.
  - `"correct"` (FHRF-C, pers. 6–10): RF pada $y_d - x_d'\hat\beta - \hat v_d$; $\hat m_d = x_d'\hat\beta + \hat f(x_d)$.
- **Algoritma** `algorithm = c("fhrf", "merf")`:
  - `"fhrf"` (default, sesuai paper): FH linear → RF (OOB untuk area tersampel, prediksi biasa untuk
    tak-sampel) → $\hat\theta_d = \hat m_d + \gamma_d(y_d - \hat m_d)$ dengan σ̂²ᵥ FH linear.
  - `"merf"` (**adaptasi sendiri, eksperimental**; DISCREPANCIES MERF-3): RF dan FH berofset diulang;
    σ̂²ᵥ diestimasi ulang dengan RF sebagai offset. Kriteria henti: perubahan σ̂²ᵥ / (σ̂²ᵥ + ψ̄) < `tol`
    (1e-3), `maxit` 25. Perubahan θ̂ tidak dipakai karena RF berbentuk tangga sehingga θ̂ terus
    berfluktuasi ±2% sd(y) (diamati pada D = 300).
- **Estimasi σ²ᵥ.** Paper memakai EM (Valdez dkk., [BELUM]). Spike S4: titik tetap EM pers. (3)–(5)
  (dengan langkah β GLS, yang tidak dicetak di paper) = ML `sae::eblupFH`. Default flexsae REML
  (sama dengan `pspline_sae`); `method = "ML"` mereproduksi paper (MERF-1).
- **Engine FH.** Pemecah profil-likelihood satu dimensi (`R/fh_lm.R`), bukan mgcv: mgcv dengan
  `s(area, bs = "re")` makan ~2,5 dtk per fit pada D = 300 (S5), sedangkan merf dan bootstrapnya
  mem-fit ulang ratusan kali. Hasilnya sama dengan `sae::eblupFH` sampai ~1e-7.
- **MSE.** Bootstrap parametrik pers. (15)–(17), fit ulang seluruh algoritma termasuk RF.
- **Asumsi/risiko.** RF dengan D kecil tidak stabil; RF memakai variabel RHS tanpa transformasi; semua
  RF dalam satu fit memakai satu seed (`ranger(seed = )`, ditarik dari stream R).
- **Batasan `pspline_sae` yang ditemukan di fase 2:** biaya mgcv O(D³) per fit (2,5 dtk pada D = 300,
  jadi bootstrap B = 200 ≈ 8 menit). Opsi: jalur cepat `.fit_fh_lm()` bila formula tanpa smooth;
  untuk smooth perlu dicari (mis. representasi marginal). Belum dikerjakan.

### 2.4 Rencana lanjutan (bukan v0.1)

Unit-level: P-spline nested error (Opsomer dkk. 2008), MERF unit-level (Krennmair & Schmid 2022;
MSE bootstrap REB Chambers & Chandra 2013 dengan varians terkoreksi Mendez & Lohr 2011), MERF
dengan auxiliary teragregasi (Krennmair, Würz & Schmid 2026), M-quantile P-spline (Chambers &
Tzavidis 2006; Pratesi, Ranalli & Salvati 2009; MSE Chambers, Chandra & Tzavidis 2011, Marchetti
dkk. 2012). Untuk metode unit-level non-linear, kovariat unit sensus wajib (ketidaksamaan Jensen).
Area-level Bayes: `bart_sae()` (Parker & Eideh 2026) setelah teks lengkapnya dibaca.

---

## 3. Pengujian terhadap SAE linear

### 3.1 Kesetaraan numerik (testthat, `skip_if_not_installed`)

| Uji | Isi | Toleransi |
|---|---|---|
| E1 | `pspline_sae()` tanpa smooth vs `sae::eblupFH` REML/ML, `emdi::fh` reml/ml, `fastsae::eblup_fh` (data `sae::milk` dan simulasi) | σ²ᵥ 1e-5; β, θ̂ 1e-6 relatif (dari S1) |
| E2 | `merf_sae()`: langkah FH berofset dengan offset tetap vs `sae::eblupFH` dengan $y_d - \text{offset}_d$ | 1e-6 |
| E3 | `merf_sae(variant = "correct")` dengan RF yang menghasilkan 0 (mis. `max.depth`/data konstan) = FH linear | 1e-6 |
| E4 | `benchmark_flexsae()` pada vektor vs `fastsae::benchmark_sae` (ratio, difference, optimal) | 1e-10 |
| E5 | vs `emdi::benchmark` (`ratio`, `raking`, `MSE_adj`; definisi dibaca dari kode emdi dulu) | 1e-10 |
| E6 | MSE terbenchmark FH-difference vs `saebenchmarking` (Steorts & Ghosh 2013) | ditentukan setelah membaca kedua sumber |

Fakta yang sudah dicek: `emdi::benchmark()` menyetel MSE ke NULL; `sae` tidak punya fungsi benchmark.

### 3.2 Simulasi (ADEMP, Morris, White & Crowther 2019 [V-meta])

- **Aims.** Akurasi titik dan interval vs FH linear; kualitas MSE; biaya komputasi; efek benchmark.
- **DGP.** $\theta_d = f(x_d) + v_d$, $y_d = \theta_d + e_d$, $e_d \sim N(0,\psi_d)$, $x_d$ baru tiap populasi atau tetap
  (dicatat). Bentuk $f$: B linear; A1 kuadratik; A2 sinus; A3 ambang; A4 interaksi $x_1x_2$; C campuran
  (2 non-linear + 1 linear + 5 noise, 2 noise berkorelasi dengan kovariat penting). **Semua bentuk saat
  ini rancangan sendiri**; diganti/dirujuk ke pengaturan simulasi Bukhari dkk. (2026) dan Parker (2024)
  setelah bagian simulasinya dibaca.
- **Variasi** (satu faktor per kali dari dasar $D = 100$): $D \in \{30, 50, 100, 300\}$; rasio
  $\sigma_v^2/\bar\psi \in \{0{,}25, 1, 4\}$; heterogenitas $\psi_d$ (homogen vs 5 kelompok); $\lambda \in \{0, 0{,}5, 1\}$;
  $v_d$ normal vs $t_3$ terskala; 10% area tak-sampel.
- **Estimand.** $\theta_d$; untuk B juga $\sigma_v^2$.
- **Methods.** `pspline_sae`, `merf_sae` (dua varian), `sae::eblupFH`, `emdi::fh`, `fastsae::eblup_fh`,
  `nlfh` (BART-FH), estimasi langsung.
- **Performance.** Bias, RMSE, RRMSE per area lalu dirangkum; coverage 95% dan lebar interval;
  skor interval (Gneiting & Raftery 2007); bias relatif dan RRMSE estimator MSE; waktu; peak RSS
  satu fit per proses (pola `bench_resources.R`).
- **Replikasi/MCSE.** Pilot $R = 100$; $R = \lceil(\mathrm{SD}/\mathrm{MCSE}_{target})^2\rceil$, target MCSE coverage ≤ 0,5
  poin persen, min 200, maks 1000; MCSE di setiap tabel; $B = 200$ untuk bootstrap MSE; evaluasi MSE
  pada subset skenario bila mahal.
- **Seed.** `RNGkind("L'Ecuyer-CMRG")`, stream per replikat, seed disimpan.
- **Hipotesis tercatat.** A: flexsae RMSE < FH linear; B: non-inferior, rasio RMSE median ≤ 1,05;
  C: lebih baik dan tidak memburuk karena noise. Tabel "skenario flexsae kalah" wajib.

### 3.3 Uji unit dan regresi

Satu file test per modul; semua kombinasi argumen; kombinasi tak didukung → `cli_abort()` di awal;
kasus tepi: area tak-sampel, $D$ sangat kecil, kovariat konstan, $\psi_d \le 0$ atau NA, σ²ᵥ = 0
(batas), tidak konvergen, MSE bootstrap dengan replikat gagal. Fixture dari `data-raw/make_fixtures.R`
menyimpan versi paket sumber. Setiap perbaikan bug disertai test yang gagal sebelum perbaikan.

---

## 4. Seleksi variabel (`select_flexsae()`)

| Model | Metode | Rujukan | Catatan |
|---|---|---|---|
| pspline | Shrinkage double-penalty (`select = TRUE`) | Marra & Wood (2011) [V-meta] | Satu fit; ambang edf heuristik, dicek di makalah |
| pspline | Stepwise IC (backward/both) dengan fit ML | Marhuenda, Morales & Pardo (2014) untuk FH; Wood, Pya & Säfken (2016) AIC terkoreksi; Vaida & Blanchard (2005) cAIC [V-meta] | Kasus linear harus sama dengan `fastsae::step_sae()` (uji regresi) |
| merf | Permutation importance pada $y_d - \hat v_d$ (OOB) | Breiman (2001); Strobl dkk. (2008) [V-meta] | Bias pada kovariat berkorelasi |
| merf | Boruta | Kursa & Rudnicki (2010) [V-meta] | `Boruta` di Suggests |

Tidak di v0.1: knockoff (Candès dkk. 2018) dan LASSO/grup (alasan di revisi sebelumnya: asumsi
pertukaran baris dan linearitas). Validasi: DGP C, ukuran proporsi replikat yang memilih semua
variabel benar, membuang semua noise, tepat himpunan benar; TPR/FPR per variabel; variasi rasio
$\sigma_v^2/\bar\psi$ (efek area besar menyamarkan sinyal kovariat). Skrip
`benchmarks/validation/select_flexsae.R`.

---

## 5. Benchmarking

> **Status 28 Okt 2026 (fase 3a selesai):** `benchmark_flexsae()` post-hoc + MSE terbenchmark bootstrap
> ada. Perubahan dari rancangan di bawah, setelah membaca kode `saebenchmarking` 0.1.0 dan `emdi` 2.2.3:
> - MSE terbenchmark **hanya bootstrap**, dihitung dengan mem-benchmark ulang tiap replikat bootstrap model
>   (replikat itu sudah fit ulang penuh; θ* dan y* kini disimpan di objek model). Rumus linearisasi §5.4
>   dan dugaan $\mathrm{Cov}(e_d,\varepsilon_g) \approx w_d\gamma_d\psi_d$ **tidak dipakai**: semua suku diestimasi dari
>   replikat. Untuk difference + target same-survey, kenaikan MSE bootstrap = g₄ Steorts & Ghosh (2013)
>   dari `saebenchmarking` (0,000078 vs 0,000077 pada milk, B = 2000; DISCREPANCIES S7).
> - Argumen: `target` numerik (tipe ii), `target` + `target_var` (tipe i eksternal), `target = "direct"`
>   (tipe i same-survey). `target_type`/`target_source` dihapus (diturunkan dari argumen).
> - **Self-benchmark (§5.3) ditunda ke fase 3b**: butirnya masih [dugaan] dan isi Wang, Fuller & Qu
>   (2008) belum dibaca penuh.

### 5.1 Notasi

Grup $g$; $w_d = N_d/\sum_{j\in g}N_j$ (`type = "mean"`) atau $w_d = N_d$ (`"total"`); galat $e = \hat\theta - \theta$,
$M = E(ee')$; target $\hat T_g$, nilai benar $T_g = \sum w_d\theta_d$, galat target $\varepsilon_g = \hat T_g - T_g$.
Metode linear: $\tilde\theta_d = \hat\theta_d + k_d(\hat T_g - w'\hat\theta)$, $\sum w_dk_d = 1$: difference, ratio, optimal
(Datta dkk. 2011; Bell, Datta & Ghosh 2013 [V-meta]); logit eksperimental (dari fastsae, tanpa rujukan).

### 5.2 `benchmark_flexsae()`

```
benchmark_flexsae(object, target, weight, group = NULL,
                  method = c("ratio", "difference", "optimal", "logit"),
                  type = c("mean", "total"),
                  target_type = c("known", "estimated"),
                  target_source = c("external", "same_survey"), target_var = NULL,
                  national_target = NULL, areas = c("all", "sampled"), bounds = NULL,
                  mse = c("linearized", "bootstrap", "none"), B = 200, seed = NULL)
```

Satu total, banyak grup, dua tahap kabupaten → provinsi → nasional (Ghosh & Steorts 2013;
Pfeffermann, Sikov & Tiller 2014 [V-meta]). Ringkasan: agregat sebelum/sesudah vs target,
penyesuaian per area, 10 area berubah terbesar, pelanggaran `bounds`, MSE sebelum/sesudah.
Area tak-sampel ikut diagregasi dan disesuaikan (`areas = "all"`, default) atau tidak
(`"sampled"`). Bobot NA/≤ 0 → `cli_abort`; bobot `mean` dinormalkan per grup dengan pesan bila
jumlahnya ≠ 1. Proporsi: peringatan bila ratio/difference keluar (0, 1); cacahan: ratio menjaga ≥ 0.

### 5.3 Self-benchmark (argumen `self_benchmark` pada fungsi model)

- `pspline_sae()`: augmentasi Wang, Fuller & Qu (2008) [V-meta, PDF dibuka] — kolom tak-berpenalti
  $z_d = w_d\psi_d$ per grup. Persamaan Henderson untuk kolom efek tetap memberi
  $\sum_d w_d(y_d - \hat\theta_d) = 0$, jadi $\sum w_d\hat\theta_d = \sum w_d y_d$ (agregat langsung area tersampel)
  [dugaan; wajib tes properti]. MSPE: You, Rao & Hidiroglou (2013) [V-meta] sebagai acuan.
- `merf_sae()`: kolom yang sama ditambahkan pada **langkah FH berofset terakhir** ($\theta_d = \hat f^{OOB}_d + z_d\alpha + v_d$);
  kendala berlaku bersyarat pada $\hat f$ [dugaan; wajib tes properti]. Jika tes gagal, merf hanya post-hoc.
- Hanya ke agregat langsung dari survei yang sama (tipe (i) same-survey), bukan ke target eksternal.

### 5.4 MSE terbenchmark

**Linearisasi (orde-1):** dari $\tilde\theta_d - \theta_d = e_d - k_dw'e + k_d\varepsilon_g$ ($k_d$ tetap):
$$\mathrm{MSE}(\tilde\theta_d) = M_{dd} - 2k_d(Mw)_d + k_d^2 w'Mw + k_d^2E\varepsilon_g^2 + 2k_d\mathrm{Cov}(e_d,\varepsilon_g) - 2k_d^2\mathrm{Cov}(w'e,\varepsilon_g).$$
Tipe (ii): tiga suku terakhir nol (tambah $k_d^2b_g^2$ bila target bias $b_g$). Tipe (i) eksternal:
$E\varepsilon_g^2$ = `target_var`, kovarians nol. Tipe (i) same-survey ($\hat T_g = \sum w_d y_d$): untuk FH dengan
parameter diketahui $\mathrm{Cov}(e_d,\varepsilon_g) \approx w_d\gamma_d\psi_d$ [dugaan]; dicocokkan dengan Steorts & Ghosh
(2013) dan You, Rao & Hidiroglou (2013) sebelum dipakai — rujukan menang bila berbeda. $M$ diestimasi
dari kovarians replikat bootstrap.

**Bootstrap ulang penuh (wajib untuk merf):** per replikat bangkitkan $\theta^*_d$, $y^*_d$; target bootstrap
$\hat T^*_g = T^*_g$ (ii), $T^*_g + \varepsilon^*$ dengan $\varepsilon^* \sim N(0,$`target_var`$)$ (i eksternal), atau $\sum w_dy^*_d$
(i same-survey); **fit ulang, prediksi, benchmark ulang**; simpan galat kuadrat terbenchmark dan tidak
terbenchmark dari replikat yang sama (perubahan MSE berpasangan, bisa naik atau turun). Dilarang
mem-benchmark kolom MSE. Keterbatasan: bootstrap di bawah model terfit tidak menangkap pengurangan
bias pada model salah spesifikasi; diukur di §5.6.

### 5.5 Output

Kelas `flexsae_benchmark`: `domain, group, weight, sampled, original, benchmarked, adjustment,
rel_adjustment_pct, target, target_type, mse, mse_benchmarked, mse_ratio, out_of_bounds`;
atribut metode, `mse_method`, `B`, seed, tabel cek agregat; `print`, `summary`, `autoplot`.

### 5.6 Pengujian benchmark

Properti $\sum_d(N_d/N)\tilde\theta_d = T$ (1e-10; logit 1e-8) untuk semua metode × `method` × `type` ×
grup/dua tahap, dan untuk self-benchmark. Akurasi bias/RMSE sebelum vs sesudah; kasus linear cocok
dengan fastsae, emdi, saebenchmarking (E4–E6). MSE: rasio MSE estimasi/empiris dan coverage untuk
target tipe (ii), (i) eksternal, (i) same-survey. Sensitivitas: target bias 1%/5%; $N_d$ salah
(lognormal 5–10%); model salah spesifikasi (FH linear pada DGP A lalu benchmark); `target_var` salah
(×0,5, ×2). Laporan wajib memuat skenario MSE naik setelah benchmark.

---

## 6. Visualisasi (`autoplot()`, pola fastsae)

`type`: `"estimates"` (langsung vs model + interval), `"scatter"`, `"mse"`, `"rse"` (ambang BPS 25/50%),
`"partial"` (pspline: kurva smooth + CI mgcv; merf: PDP Friedman 2001, ICE Goldstein dkk. 2015, ALE Apley
& Zhu 2020), `"importance"`, `"residuals"` (residual terstandar area, QQ $\hat v_d$), `"caterpillar"`
($\hat v_d$ ± interval dari varians kondisional atau replikat bootstrap), `"map"` (`sf`, Suggests);
`autoplot.list()` untuk perbandingan metode; `autoplot.flexsae_benchmark()`.

---

## 7. Dokumentasi pkgdown

`_pkgdown.yml` pola fastsae (bootstrap 5, flatly). Vignette: getting-started; area-level-nonlinear
(model, MSE, kapan memakai); variable-selection; benchmarking; linear-comparison (E1–E6 + ADEMP,
termasuk skenario kalah); visualization. Angka dari `inst/extdata/*.rds` lewat
`benchmarks/update_benchmark_docs.R`.

---

## 8. Arsitektur

```
R/  flexsae-package.R utils.R checks.R
    pspline_sae.R merf_sae.R fh_gam.R (FH via mgcv, pspline) fh_lm.R (FH satu komponen varians, merf)
    mse_bootstrap.R benchmark.R self_benchmark.R
    select.R importance.R partial.R methods.R autoplot.R compare.R sim_data.R data.R
tests/testthat/  test_<modul>.R, helper-*.R, fixtures/
data-raw/  benchmarks/  technical-report/  vignettes/
```

Objek: kelas `c("flexsae_pspline" | "flexsae_merf", "flexsae")`; `df_est` (`domain, direct, vardir, est,
synthetic, random_effect, gamma, mse, rse, ci_lower, ci_upper, sampled`), `fit`, `sigma2_v`, `formula`, `call`, `data`, `model`,
`mse_method`, `B`, `seed`, `convergence`, `iterations`, `self_benchmark`. Kolom yang tidak berlaku = NA.

DESCRIPTION target: Imports `cli, ggplot2, lifecycle, mgcv, ranger, rlang, stats, utils`; Suggests
`sae, emdi, fastsae, saebenchmarking, nlfh, Boruta, sf, ps, withr, testthat (>= 3.0.0), knitr,
rmarkdown, covr`. GPL (>= 3).

---

## 9. Tahapan

| Fase | Isi | Selesai bila |
|---|---|---|
| 0 ✅ | Repo publik, kerangka package, CI; spike S1 lulus | R CMD check bersih pada kerangka |
| 1 ✅ | inti FH mgcv (`R/fh_gam.R`) + `pspline_sae()` (titik, bootstrap MSE, S3), E1 vs sae/fastsae/emdi | test hijau, check 0 ERROR/WARNING |
| 2 ✅ | `merf_sae()` dua varian × dua algoritma + bootstrap MSE, E2–E3, pemecah FH cepat | idem |
| 3a ✅ | `benchmark_flexsae()` post-hoc, MSE terbenchmark bootstrap, E4–E6, tes properti | idem |
| 3b | self-benchmark (WFQ) untuk `pspline_sae`/`merf_sae` setelah paper dibaca | idem |
| 4 | `select_flexsae()`, importance, partial, `autoplot` | idem + validasi seleksi pilot |
| 5 | Simulasi ADEMP (paket terinstal), vignette, pkgdown | tabel dengan MCSE, situs terbangun |
| 6 | Rencana lanjutan §2.4 | – |

Setiap fase: `devtools::document()` → `devtools::test()` → `R CMD check --as-cran`; hasil dilaporkan
apa adanya. win-builder/CRAN/tag hanya setelah konfirmasi.

---

## 10. Keputusan

Sudah diputuskan (28 Okt 2026): area-level dulu; repo publik; model `*_sae`, utilitas `*_flexsae`;
v0.1 = `pspline_sae` + `merf_sae`; maintainer seperti fastsae.

Masih terbuka: sisa prompt yang terpotong; fastsae sebagai Suggests (usulan) atau Imports; DGP
rancangan sendiri vs menunggu pengaturan simulasi Bukhari dkk./Parker dibaca.

## 11. Risiko

$D$ kecil (overfit RF dan spline); rujukan FH P-spline frekuentis belum terverifikasi; venue Bukhari
dkk.; derivasi sendiri (self-benchmark merf, kovarians same-survey) bisa salah — semuanya diberi tes;
biaya bootstrap merf × replikasi simulasi.
