# One-command test36 pipeline

AA-MDの保存済み軌跡から、CGデータ作成、学習、構造生成、近似CG-MD、RDF解析までを一括実行する。

```bash
DM2_ROOT="$PWD/DM2" DEVICE=mps bash run_test36_all.sh
```

主な出力は`test36-aa/long-run/{dataset,train,generated,cgmd,cgmd-analysis}`。
既存の完了ファイルを再利用し、途中失敗後も同じコマンドで続行できる。ただし壊れた途中出力は新しいディレクトリへ移してから再実行する。

設定例：

```bash
DEVICE=cuda TRAIN_TIME_HOURS=11.5 UPDATES=30000 \
  DATASET_DIR=/work/project/test36-data \
  TRAIN_DIR=/work/project/test36-train \
  bash run_test36_all.sh
```

`test36.py`のモデルと拡散アルゴリズムはtest32と同じで、学習条件をモデル入力にはしていない。出力は厳密な自由エネルギー力場ではなく、`sigma-ref`に基づく近似CG-MDである。
