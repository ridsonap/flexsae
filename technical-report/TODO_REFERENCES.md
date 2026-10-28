# TODO_REFERENCES — rujukan yang BELUM terverifikasi atau perlu catatan

Aturan (sama dengan fastsae): tidak ada DOI/entri yang dikarang. Yang tidak dapat
dikonfirmasi ke Crossref, situs penerbit, atau PDF resmi dicatat di sini dan **tidak
dikutip** di kode, Rd, vignette, maupun laporan.

Tanggal verifikasi: 28 Okt 2026. Entri terverifikasi ada di `references.bib` (68 entri).
Lingkup v0.1 kini **area-level** saja (lihat `PLAN.md`); rujukan unit-level tetap disimpan
untuk rencana lanjutan.
Cara verifikasi: `api.crossref.org/works/<doi>` untuk entri ber-DOI; PDF di
`www150.statcan.gc.ca` (katalog 12-001-X) untuk *Survey Methodology*; situs
*Statistica Sinica* untuk volume/halaman Steorts & Ghosh (2013).

## 1. Belum terverifikasi → tidak boleh dikutip

| Rujukan | Status | Alasan / langkah berikut |
|---|---|---|
| Giusti, Marchetti, Pratesi & Salvati (2012), "Semiparametric Fay-Herriot model using penalized splines", *J. Indian Soc. Agric. Statist.* 66(1):1–14 | **BELUM** | Hanya ditemukan sebagai sitasi sekunder (daftar pustaka Parker 2024 dan Parker & Eideh 2026, ResearchGate, Google Scholar). Tidak ada DOI di Crossref; halaman penerbit belum dibuka. Ini satu-satunya rujukan *frekuentis* untuk FH P-spline area-level, jadi perlu diverifikasi sebelum `fsae_fh_pspline()` dirilis. Langkah: cari di situs ISAS (isas.org.in) atau minta PDF ke penulis. |
| Isaki, Tsay & Fuller (2000), "Estimation of census adjustment factors", *Survey Methodology* 26:31–42 | **BELUM** | Hanya sitasi sekunder (daftar pustaka Pfeffermann, Sikov & Tiller 2014). PDF StatCan belum dibuka. |
| You, Rao & Dick (2004), *Statistics in Transition* 6:631–640 (benchmarking HB) | **BELUM** | Hanya sitasi sekunder (Ugarte dkk. 2009; Pfeffermann dkk. 2014). |
| Messer & Schmid (2024), gradient boosting untuk data hierarkis di SAE, arXiv:2406.04256 | **BELUM** | Hanya sitasi sekunder di Parker & Eideh (2026), yang menyebutnya "preprint: not peer reviewed". Tidak dipakai. |
| Fisher, Rudin & Dominici (2019), model reliance / permutation importance, *JMLR* 20 | **BELUM** | Tidak dicek. Permutation importance cukup dirujuk ke Breiman (2001) dan Strobl dkk. (2008) yang sudah terverifikasi. |
| Hajjem, Bellavance & Larocque (2011), mixed effects regression trees | **BELUM** | Tidak dicek; yang dipakai adalah Hajjem dkk. (2014, MERF), terverifikasi. |
| Valdez dkk., algoritma EM khusus Fay-Herriot (rujukan [22] di Bukhari dkk. 2026) | **BELUM** | Hanya diketahui dari teks Bukhari dkk.; penulis lengkap, judul, dan tahun belum dicari. `merf_sae()` memakai REML FH dengan offset (Fay & Herriot 1979; spike S1), bukan EM ini. |
| Das, Jiang & Rao (2004), MSE EBLUP untuk model campuran umum | **BELUM** | Disebut di PLAN sebagai *kemungkinan* dasar MSE analitik multi-komponen varians; belum dicari. Sampai diverifikasi, MSE FH-spline memakai bootstrap parametrik. |

## 2. Terverifikasi metadatanya, tetapi ISI belum dibaca penuh

Boleh dikutip sebagai rujukan umum. Persamaan atau prosedur spesifik **tidak boleh**
diimplementasikan sebelum teks lengkapnya dibaca dan nomor persamaannya dicatat
(aturan skill `r-sae-package` §1).

| Entri | Yang sudah dibaca | Yang belum |
|---|---|---|
| `bukhari2026random` | Halaman penerbit (iieta.org): abstrak, persamaan (1)–(17): varian FHRF-C dan FHRF-R, EM untuk σ²ᵥ, bootstrap parametrik (15)–(17) | Apakah RF memakai prediksi OOB atau in-sample; pengaturan RF; detail DGP simulasi. **Catatan venue:** *Ingénierie des Systèmes d'Information* (IIETA) berprofil rendah untuk metodologi statistik; dasar metodologis MERF tetap `hajjem2014mixed` dan `krennmair2022flexible`. |
| `bukhari2025modified` | Metadata (prosiding *BIO Web of Conferences*) | Isi |
| `opsomer2008nonparametric` | Abstrak (Crossref): ada hasil teoretis MSE prediksi dan bootstrap non-parametrik | Rumus MSE analitik dan langkah bootstrap |
| `krennmair2022flexible` | Versi arXiv 2201.10933v2: model, algoritma EM, prediktor, bootstrap REB, skenario simulasi | Versi terbit JRSS C (cek apakah ada perubahan dari preprint) |
| `parker2026bartfh` | Metadata, daftar pustaka | Model, prior, dan cara inferensi (halaman OUP tidak memuat teks lengkap) |
| `parker2024nonlinear` | Persamaan model (PDF NSF-PAR) | Prior lengkap dan pemilihan jumlah neuron |
| `sugasawa2019hierarchical`, `polettini2017generalised` | Judul/metadata | Isi |
| `ugarte2009benchmarked` | Abstrak: estimator terestriksi + MSE prediksi diturunkan | Apakah model unit-level; bentuk MSE |
| `steorts2013estimation` | Abstrak: FH, kenaikan MSE akibat benchmark O(m⁻¹), estimator MSE orde-2 dan bootstrap parametrik | Persamaan |
| `ghosh2015benchmarked`, `ghosh2013two`, `bell2013benchmarking`, `datta2011bayesian` | Metadata | Persamaan |
| `wang2008small` | PDF: abstrak, model augmentasi | Rumus MSE prediktor terkalibrasi |
| `you2013performance`, `stefan2021small` | PDF: abstrak | Rumus MSPE |
| `chambers2011bias`, `marchetti2012nonparametric` | Abstrak | Rumus |

## 3. Koreksi terhadap prompt dan terhadap ingatan saya

Rujukan yang disebut di prompt (dari ingatan) **semuanya ternyata ada**; detailnya:

- Hajjem dkk.: terbit daring 2012, volume cetak **2014** (84(6):1313–1328). Dipakai 2014.
- Sela & Simonoff: daring 2011, cetak **2012** (86(2):169–207). Dipakai 2012.
- Datta, Ghosh, Steorts & Maples: daring 2010, cetak **2011** (TEST 20(3):574–588).
  Dipakai 2011. (Catatan: `fastsae/technical-report/references.bib` memuat entri ganda
  `datta2010bayesian` dan `datta2011bayesian` dengan DOI yang sama.)
- "Pratesi dkk." → Pratesi, Ranalli & Salvati (2009), *J. Nonparametric Statist.* 21(3).
- "Wood" → Wood (2011, JRSS B; daring 2010) dan Wood (2017, buku edisi ke-2).
- You & Rao (2002) **punya DOI** `10.2307/3316146` (Crossref). `fastsae/technical-report/TODO_REFERENCES.md`
  §5 menyebutnya "tanpa DOI terdaftar"; itu keliru, meskipun bib fastsae sudah memuat DOI-nya.

DOI dari ingatan saya yang **salah** (diperbaiki lewat pencarian Crossref):

| Rujukan | DOI salah (milik makalah lain) | DOI benar |
|---|---|---|
| Lombardía, López-Vizcaíno & Rueda (2017) | 10.1111/rssa.12309 | 10.1111/rssa.12300 |
| Pfeffermann & Barnard (1991) | 10.1080/07350015.1991.10509830 | 10.1080/07350015.1991.10509828 |
| Sperlich & Lombardía (2010) | 10.1080/10485250903459718 (tidak ada) | 10.1080/10485250903311607 |
| Rao, Sinha & Dumitrescu (2014) | 10.1002/cjs.11204 | 10.1002/cjs.11199 |

## 4. Metadata terverifikasi dengan catatan

- `krennmair2026analysing`: Crossref "issued" 2025, volume cetak 75(1) tahun 2026. Dipakai 2026.
- `parker2026bartfh`: *advance article* JSSAM, belum ada volume/halaman. Kolom `pages` berisi
  nomor artikel `smaf050`. Perbarui setelah terbit cetak.
- `molina2015sae`: Crossref hanya mencatat halaman awal (81); rentang 81–98 dari situs R Journal.
- `steorts2013estimation`: Crossref tanpa volume/halaman; 23(2):749–767 dari situs Statistica Sinica.
- Entri *Survey Methodology* (`wang2008small`, `you2013performance`, `stefan2021small`,
  `chambers2011bias`) tidak punya DOI; dicek di PDF StatCan.
- `butar2003measures`: nama penulis mengikuti Crossref (`Butar, Ferry Butar`).

## 5. Paket R terkait (bukan rujukan bibliografis, dicek di CRAN 28 Okt 2026)

| Paket | Status | Relevansi |
|---|---|---|
| `flexsae` | **Belum ada** di CRAN (halaman 404, tidak ada arsip) | Nama tersedia |
| `SAEforest` | **Diarsipkan** 2023-08-10 ("issues were not corrected in time"), versi terakhir 1.0.0 | Hanya untuk skrip `benchmarks/` dan fixture; tidak boleh di Suggests |
| `nlfh` 0.1.0 (P. Parker, 2026-06-12) | Di CRAN, berbasis `dbarts` | Pembanding BART-FH area-level |
| `saebenchmarking` (2026-09-28) | Di CRAN | Pembanding MSE terbenchmark FH (Steorts & Ghosh 2013) |
| `saekernel` 0.1.1 | Di CRAN (email maintainer tidak terkirim) | FH Nadaraya–Watson; rujukannya tidak diverifikasi, tidak dipakai |
| `sae` 1.3, `emdi` 2.2.3, `mgcv` 1.9-4, `ranger` 0.18.0, `dbarts` 0.9-33 | Terpasang lokal | Pembanding/engine |
