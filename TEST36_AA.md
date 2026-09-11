# test36用：粘土＋水＋Na/Caの全原子MD

`test36_aa.py` は実際の粘土単位格子と古典力場から、明示的な水と
Na/Caを含む全原子MDを実行する。`test32.py`/`test36.py` の拡散モデルは
変更しない。乱数による水の配置は**初期配置だけ**で、学習に渡す軌跡は
CLAYFF系の力とOpenMM積分器で進めた座標。

## 今回の条件

| 項目 | 初期値／設定 |
|---|---|
| 粘土 | 理想化cis型モンモリロナイト、Si₈Al₃.₅Mg₀.₅O₂₀(OH)₄ |
| セル | 4×4単位格子×2層、三次元周期境界 |
| 粘土原子 | 1,280原子（Si256、Al112、Mg16、O768、H128） |
| 水 | rigid SPC 384分子＝1,152原子 |
| 対イオン | Na⁺8個＋Ca²⁺4個、Na/Caの電荷当量1:1 |
| 全原子数 | 2,444 |
| 総電荷 | 粘土−16e＋対イオン16e＝0（丸め誤差~4e−13e） |
| ボックス | 20.64×35.864×31.0 Å |
| 初期層間隔 | 15.5 Å（層厚も含むd-spacing） |
| 条件 | 300 K、NVT、追加塩なし、全粘土原子可動 |
| 積分 | Langevin Middle、0.5 fs、摩擦1 ps⁻¹ |
| 静電相互作用 | PME、許容誤差1e−5、実空間カットオフ10 Å |
| LJ | Lorentz–Berthelot混合、10 Å打ち切り、tail補正・switchなし |
| 試走 | 最小化→昇温2 ps→緩和5 ps→保存対象20 ps |
| 保存間隔 | 200 step＝0.1 ps、計200フレーム |

粘土格子は無限周期の2枚のシートで、端面や有限サイズのtactoidはない。
Mg置換は規則的な配置で、実試料の化学組成・ランダム置換を再現した系ではない。
水384分子・層間隔15.5 Åは**選んだ初期条件**で、膨潤平衡から決まった値ではない。
ボックスは固定だが、各シートの位置・形状・個々の層間距離は自由に動く。

## 実行結果の場所

ローカルCPU（4 threads）で27 psが正常終了。最小化を含むMD実行時間は
258.05秒（約4分18秒）、保存対象20 psの平均温度は300.23 K、
フレーム間標準偏差は5.03 K。最終フレームでSiの4配位、Al/Mgの6配位は
それぞれ100%（解析カットオフは2.2/2.6/2.7 Å）。
CG入力は200フレーム×1,292サイトで、先頭・末尾フレームのAA→CG対応と
データハッシュを検査済み。学習はまだ実行していない。

- `test36-aa/initial-v2/`: 本試走の初期構造、`system.top`、`mapping.json`、出典。
- `test36-aa/pilot/production.extxyz`: 全原子の保存対象軌跡（昇温・緩和を除く）。
- `test36-aa/pilot/thermo.csv`: 全段階の温度・エネルギー・計算速度。
- `test36-aa/pilot/run_metrics.json`: 完了状況と試走平均温度など。
- `test36-aa/pilot/analysis.json`: 配位数、水の拘束誤差、温度・エネルギーのブロック平均。
- `test36-aa/pilot/latest.chk`: OpenMMのバイナリチェックポイント。
- `test36-aa/pilot/latest.state.xml`: 可搬な座標・速度などの状態。
- `test36-aa/dataset-pilot/`: `test36.py prepare` が作ったCG入力。
- `test36-aa/pilot.log`: tmux内のMD、変換、解析のログ。

`initial/` と `smoke-cpu/` は最初の0.1 psのソフトウェア試走用。
水分子を周期境界で分断しない処理を追加した後の `initial-v2/` を使う。
最初のsmoke軌跡は学習用データに混ぜない。

## 出典と変更点

[ClayCode公開リポジトリ](https://github.com/Erastova-group/ClayCode) の
commit `fb52753bbad9c21369442589ac66b16a3bdb24e2` を固定。
元データと付属MITライセンスは `test36-aa/sources/`、検査用SHA-256は
`test36-aa/source_manifest.json` に保存。ビルダーはハッシュが違えば停止する。

- CD21のC2000（未置換）、C2013/C2014（八面体Al→Mg置換）を組み合わせる。
  元の各原子の微小な電荷補正も維持する。Mgの6配位酸素と
  `obos`/`ohs`の一致をセル境界も含めて検査する。
- 粘土は配布の `ClayFF_Fe.ff`（今回はFeなし）、水は配布の `spc.itp`、
  Na/Caは配布の `Ions.ff` のIODパラメータを使う。
- 元のSiタイプ `st` に原子番号16と記載されているため、**元素ラベルだけ14に修正**。
  質量・電荷・LJ・OH結合パラメータは変えない。
- 水はO–H=1 Å、H–H=1.633 Åで剛体。粘土OHは元の調和結合で可動。
  分子内除外は水の3組と粘土のOH結合のみ。粘土のSi–O、Al–Oなどは非結合相互作用。
- ClayCode本体の構築・GROMACS実行コードを動かしたわけではない。
  **配布データを使う独自ビルダー＋OpenMM**であり、論文のMDプロトコル完全再現ではない。

参考文献：

1. [Pollak et al., Modeling Realistic Clay Systems with ClayCode (2024)](https://doi.org/10.1021/acs.jctc.4c00987)
2. [Cygan et al., Advances in Clayff Molecular Simulation of Layered and Nanoporous Materials (2021)](https://doi.org/10.1021/acs.jpcc.1c04600)
3. [Smith et al., IOD等のイオンパラメータ (2023)](https://doi.org/10.1021/acs.jctc.2c01255)
4. [OpenMM: GROMACS入力・MD実行](https://docs.openmm.org/latest/userguide/application/02_running_sims.html)

## 再実行

依存：Python、OpenMM 8.2、NumPy、SciPy、ASE。CG変換は別途 `test36.py`
のPyTorch/graphite依存が必要。ローカルの既存Pixi環境で動作確認。

```bash
../.pixi/envs/default/bin/python test36_aa.py build \
  --output test36-aa/initial-new

../.pixi/envs/default/bin/python test36_aa.py run \
  --input test36-aa/initial-new --output test36-aa/run-new \
  --platform CPU --threads 4 --warmup-ps 2 --equilibration-ps 5 --production-ps 20

../.pixi/envs/default/bin/python test36.py prepare \
  --trajectory test36-aa/run-new/production.extxyz \
  --mapping test36-aa/run-new/mapping.json --output test36-aa/dataset-new

../.pixi/envs/default/bin/python test36_aa.py analyze \
  --input test36-aa/initial-new --run test36-aa/run-new
```

出力先は新規または空のディレクトリに限る。既存結果を上書きしない。
一括実行は `bash run_test36_aa_local.sh`。`AA_INPUT`, `AA_OUTPUT`,
`CG_DATASET`, `PRODUCTION_PS`, `EQUILIBRATION_PS`, `AA_THREADS` 等を変更可能。

```bash
tmux new-session -d -s test36-aa-new \
  'AA_OUTPUT=test36-aa/run-new CG_DATASET=test36-aa/dataset-new bash run_test36_aa_local.sh > test36-aa/run-new.log 2>&1'
```

## test36へのマッピング

今回は**水だけを消去する最小限の溶媒粗視化**。粘土のSi/Al/Mg/O/OHのHと
Na/Caはすべて残す。全2,444原子→1,292サイト、10種類の力場タイプ。
粒子の重心まとめ・Hの質量移動・剛体板への変換はしていない。
表面Oと水Oはタイプと固定原子IDで区別し、酸素を一括削除しない。

`positions.npy` は `(200,1292,3)`、`cells.npy` は `(200,3,3)`、単位はÅ。
`metadata.json` に条件とマッピングとハッシュを記録する。
これで `test36.py train --dataset test36-aa/dataset-pilot ...` に渡せるが、
**短い同一軌跡からの相関したデータであり、検証済み平衡学習データではない**。
全生成物で `source_is_equilibrium=false` とする。

## スパコン

全原子MDだけなら `test36_aa.py` と `initial-v2/` をコピーすればよく、
実行にはDM2/PyTorchも元の単位格子集も不要。初期構造から作り直す場合は
`source_manifest.json` と元データ／ライセンスもコピーする。
OpenMMでCUDA対応の環境を用意して `--platform CUDA` を選ぶ。
`run_test36_aa.pbs` は環境のPythonを使う雛形で、キュー名・資源・モジュールは
導入先に合わせる。CUDA/PBS実行はこのMacでは未検証。

```bash
qsub -P <ProjectGroup_ID> \
  -v AA_INPUT=/path/to/initial-v2,AA_OUTPUT=/path/to/aa-long,AA_PYTHON=/path/to/env/bin/python,PRODUCTION_PS=1000,EQUILIBRATION_PS=500 \
  run_test36_aa.pbs
```

例の500 ps緩和＋1 ns保存は「必ず平衡に達する時間」ではない。
時系列の収束・独立な初期条件・水/イオン分布・セルサイズ/カットオフ依存性を
別途確認する。NPTによる膨潤平衡、塩濃度系列、剪断、圧力検証は未実装。

`latest.chk` は同一OpenMM版・プラットフォーム等に依存する。
`latest.state.xml` は移植しやすいが乱数状態の厳密な継続にはならない。
現CLIは途中再開を実装していないので、未完了を同じ出力に再実行しない。
SIGTERM/時間制限は保存間隔の境界で状態を保存してexit 75。
最小化中は割込み/時間予算のチェックがMDループに入るまで遅れる。

## テスト

```bash
../.pixi/envs/default/bin/python -m unittest test36_aa_checks -v
```

組成・水分子形状・中性化・マッピング・ハッシュ・力場の単位・除外・
周期並進不変性・短いMD・既存出力保護を検査する。これらは実験への一致や
長時間の平衡性を証明するテストではない。
