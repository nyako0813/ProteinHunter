# Claude Code 実装指示書: domain_family_map v3(既存カテゴリの過剰一致監査・修正)

## 背景・目的

domain_family_map v2実装後、MA_0795・MA_0826(radical SAM/Elp3型)をクエリに実行したところ、`radical_sam x electron_carrier`ルールがMA_4170・MA_1812という候補で発火した。ここでユーザーから「ルールセットが少ないためにelectron_carrierが上位に来ているだけではないか。それを検証するには、まず既存カテゴリの定義精度を監査すべきではないか」という指摘があった。

調査の結果、**この指摘は正しく、しかも具体的なバグが見つかった**。`atp_dependent_activator`で以前見つけた「HUPスーパーファミリーのような広い構造フォールドのInterPro/SUPFAMアクセッションを判定に使うと、無関係な酵素まで拾ってしまう」という問題(`claude_code_instructions_domain_family_map_v2.md`の問題2)と**全く同じ種類のバグ**が、`electron_carrier`と`cysteine_desulfurase`にも存在していた。

`match_fields`機構(v2で新設済み)を使えば、コード変更なしで既存カテゴリの定義(YAML)を修正するだけで直せる。今回のタスクは**新規カテゴリの追加ではなく、既存10カテゴリ全件の定義精度監査と、見つかった過剰一致の修正**。

## 監査方法

UniProtバルクエクスポート(`uniprotkb_methanosarcina_acetivorans_2026_09_06.json`)に対し、各カテゴリについて「pfamアクセッションのみで一致する候補数」と「pfam+interpro+supfamすべてで一致する候補数(=現状の`categories_for()`の挙動)」を比較した。差が大きいカテゴリほど、フォールドレベルの広いInterPro/SUPFAMアクセッションに引きずられて無関係な候補まで拾っている可能性が高い。

## 監査結果

| カテゴリ | pfamのみ | pfam+interpro+supfam(現状) | 差分 | 判定 |
|---|---|---|---|---|
| `electron_carrier` | 3 | 32 | +29 | **要修正**(重大) |
| `atp_dependent_activator` | 1 | 40 | +39 | 修正済み(v2で`match_fields: [pfam]`適用済み、対応不要) |
| `sulfur_carrier` | 5 | 5 | 0 | 問題なし |
| `cysteine_desulfurase` | 6 | 23 | +17 | **要修正**(重大) |
| `rhodanese_sulfurtransferase` | 5 | 6 | +1 | 要修正(軽微) |
| `thif_moeb_activator` | 1 | 1 | 0 | 問題なし |
| `tusa_sulfur_relay` | 1 | 1 | 0 | 問題なし |
| `dsre_sulfur_relay` | 3 | 6 | +3 | 要修正(pfamリストへの追加が適切) |
| `radical_sam` | 40 | 40 | 0 | 問題なし |
| `iron_sulfur_cluster` | 19 | 19 | 0 | 問題なし(`exclusive: true`のためinterpro/supfamは元々未使用) |

## 修正詳細

### 1. `electron_carrier`(重大・修正必須)

`IPR029039`(Flavoprotein-like_sf)と`SSF52218`(Flavoproteins)は、フラビン結合フォールドを共有する非常に広い構造スーパーファミリーを指す。実データでは、本来の狙いである狭義のフラボドキシン(`PF12724`、MA_0361/MA_0363/MA_2078の3件)に加えて、以下の29件が誤って拾われている:

- `PF03358`(FMN_red、NADPH依存FMN reductase様/鉄硫黄フラボプロテイン)系: 21件(MA_1812, MA_4356, MA_3966, MA_3740, MA_3658等)
- `PF02525`(quinone reductase)系: 2件
- `PF00258`(canonical flavodoxin、これも実は妥当な可能性があるが今回は狭義に倣いpfamリストへの追加はしない。対象外セクション参照): 2件
- `PF12682`, `PF19583`系: 3件
- Pfamヒットが空でinterpro/supfamのみ一致(MA_4170含む): 1件

MA_0795/MA_0826クエリで上位に来た`MA_1812`(PF03358)・`MA_4170`(Pfamヒットなし)は、この過剰一致によって`electron_carrier`に誤分類されていた可能性が高い。

**修正**: `match_fields: [pfam]`を追加。`pfam: [PF12724]`は変更しない。interpro/supfamのアクセッションはYAMLに参考情報として残す。

```yaml
electron_carrier:
  pfam: [PF12724]
  interpro: [IPR026816, IPR029039]   # 参考情報のみ、判定には使わない(2026-09-09: match_fields追加により無効化)
  supfam: [SSF52218]                  # 同上
  match_fields: [pfam]
  note: "flavodoxin/flavoprotein様の電子伝達体。MA_0361/MA_0363で確認済み(2026-09-09)。interpro/supfamはFlavoprotein-like_sfという広い構造スーパーファミリーを指し無関係な酵素(NADPH依存FMN reductase等)まで拾うため、match_fields: [pfam]でPF12724のみに限定(2026-09-09監査で発見・修正)。"
```

修正後、`electron_carrier`は3候補(MA_0361, MA_0363, MA_2078)のみに一致する。

### 2. `cysteine_desulfurase`(重大・修正必須)

`IPR010970`/`IPR016454`/`SSF53383`(PLP-dependent transferases)は、ビタミンB6依存性酵素全般という非常に広いスーパーファミリーを指す。実データでは、狙いのシステイン脱硫酵素(`PF00266`、6件)に加えて、以下の17件が誤って拾われている:

- アミノトランスフェラーゼ類(aspB, tyrB, hisC, argD等、`PF00155`/`PF00202`系): 多数
- デカルボキシラーゼ類(グルタミン酸デカルボキシラーゼ等、`PF00282`): 2件
- セリンヒドロキシメチルトランスフェラーゼ(glyA, `PF00464`)、メチオニンガンマリアーゼ(mgl, `PF01053`)等

**修正**: `match_fields: [pfam]`を追加。`pfam: [PF00266]`は変更しない。

```yaml
cysteine_desulfurase:
  pfam: [PF00266]   # Aminotran_5
  interpro: [IPR010240, IPR010970, IPR016454]   # 参考情報のみ、判定には使わない(2026-09-09: match_fields追加により無効化)
  supfam: [SSF53383]                              # 同上
  match_fields: [pfam]
  note: "システイン脱硫酵素(IscS/SufS型)。実データ: MA_3264(iscS), MA_2718(iscS), MA_0236, MA_0808。interpro/supfamはPLP依存酵素全般という広いスーパーファミリーを指し、無関係なアミノトランスフェラーゼ・デカルボキシラーゼまで拾うため、match_fields: [pfam]でPF00266のみに限定(2026-09-09監査で発見・修正)。"
```

**既知の残存課題**: `PF00266`(Aminotran_5)自体が、システイン脱硫酵素(IscS/SufS型)と一部のクラスVアミノトランスフェラーゼ(`MA_1950`, `MA_1816`/aspC)の両方を含むやや広いPfamファミリーであるため、pfamのみに限定してもなお2件程度の偽陽性が残る可能性がある。これはPfamファミリー自体の粒度の限界であり、`match_fields`では解決できない。将来的に必要であれば、配列レベルのモチーフ判定等、別の仕組みでの絞り込みを検討する(今回は対象外)。

### 3. `rhodanese_sulfurtransferase`(軽微・修正推奨)

`IPR001763`/`SSF52821`により、`MA_4158`(A-type ATP synthase subunit A、硫黄転移とは無関係)が1件誤って拾われている。

**修正**: `match_fields: [pfam]`を追加。`pfam: [PF00581]`は変更しない。

```yaml
rhodanese_sulfurtransferase:
  pfam: [PF00581]
  interpro: [IPR001763]   # 参考情報のみ、判定には使わない(2026-09-09: match_fields追加により無効化)
  supfam: [SSF52821]       # 同上
  match_fields: [pfam]
  note: "ロダネーゼ様硫黄転移ドメイン。実データ: MA_3178, MA_2500, MA_0746。interpro/supfamがA-type ATP synthase subunit A(MA_4158)を誤って拾うため、match_fields: [pfam]で限定(2026-09-09監査で発見・修正)。"
```

### 4. `dsre_sulfur_relay`(要修正・pfamリストへの追加)

`IPR027396`/`IPR003787`により、狙いの3件(MA_3926, MA_3927, MA_3072、いずれも`PF02635`)に加えて3件が拾われるが、内訳は以下の通り**一様ではない**:

- `MA_3928`(DsrH like protein、`PF04077`)は**真の陽性**である可能性が高い。MA_3926・MA_3927・MA_3928は遺伝子番号が3つ連続しており(オペロンの可能性)、大腸菌のTusBCD(TusB+TusC+TusD、DsrEFH命名ではDsrE+DsrF+DsrH)ヘテロ六量体複合体構造(Structure誌2006、既出)と同じ3遺伝子構成に対応する。これはPfamファミリーの粒度の問題ではなく、**関連するが別のPfamファミリー(DsrH, PF04077)を`pfam`リストに追加すべきケース**。
- `MA_0809`(NADH dehydrogenase、`PF13686`)・`MA_0398`(Uncharacterized protein、`PF13686`)は無関係と判断し除外する。

**修正**: `match_fields: [pfam]`を追加した上で、`pfam`リストに`PF04077`(DsrH)を追加する。

```yaml
dsre_sulfur_relay:
  pfam: [PF02635, PF04077]   # DsrE, DsrH(2026-09-09追加: MA_3928がMA_3926/3927と連続する3遺伝子オペロン構成のため)
  interpro: [IPR027396, IPR003787, IPR007215]   # 参考情報のみ、判定には使わない(2026-09-09: match_fields追加により無効化)
  match_fields: [pfam]
  note: >
    TusBCD/DsrEFH型硫黄リレータンパク質(TusAから硫黄を受け取る側)。
    実データ: MA_3926, MA_3927, MA_3928(3遺伝子連続、TusBCD型オペロンの可能性)、MA_3072。いずれも現状「Uncharacterized protein」/「DsrH like protein」、機能未確定。
    interpro/supfamがNADH dehydrogenase等の無関係タンパク質(MA_0809, MA_0398)を拾うため、match_fields: [pfam]で限定した上、
    真の陽性であるMA_3928(DsrH, PF04077)をpfamリストに明示的に追加(2026-09-09監査で発見・修正)。
```

### 5. 修正不要なカテゴリ

- `atp_dependent_activator`: v2で`match_fields: [pfam]`を既に適用済み。対応不要。
- `sulfur_carrier`, `thif_moeb_activator`, `tusa_sulfur_relay`: pfamのみとpfam+interpro+supfamで一致数に差がなく、問題なし。
- `radical_sam`, `iron_sulfur_cluster`: interpro/supfamを元々設定していない(`iron_sulfur_cluster`は`exclusive: true`のため設計上pfamのみで判定)。問題なし。

## 実装タスク

### 1. `config/domain_family_map.v1.yaml`の更新

上記の`electron_carrier`・`cysteine_desulfurase`・`rhodanese_sulfurtransferase`・`dsre_sulfur_relay`の4カテゴリを、上記のYAMLブロックの内容に差し替える。他の6カテゴリ(`atp_dependent_activator`, `sulfur_carrier`, `thif_moeb_activator`, `tusa_sulfur_relay`, `radical_sam`, `iron_sulfur_cluster`)は変更不要。`category_rules`も変更不要(カテゴリ名は変わらない)。

コード変更は不要(`match_fields`フィールドはv2で既に実装済み)。

### 2. テスト

`tests/test_domain_family_map.py`の`test_shipped_v1_domain_family_map_loads`に、以下のアサーションを追加:

```python
assert loaded.categories["electron_carrier"].match_fields == frozenset({"pfam"})
assert loaded.categories["cysteine_desulfurase"].match_fields == frozenset({"pfam"})
assert loaded.categories["rhodanese_sulfurtransferase"].match_fields == frozenset({"pfam"})
assert loaded.categories["dsre_sulfur_relay"].match_fields == frozenset({"pfam"})
assert loaded.categories["dsre_sulfur_relay"].pfam == frozenset({"PF02635", "PF04077"})
```

`tests/test_interaction_scoring.py`に、以下の非回帰テストを追加することを推奨(必須ではない):

- `electron_carrier`カテゴリが、PF03358(NADPH依存FMN reductase様)のみを持つ候補には一致せず、PF12724(真のフラボドキシン)を持つ候補には一致することを確認するテスト。

### 3. 動作確認

```bash
cd ~/projects/ProteinHunter
python -m pytest -q
```

その後、MA_0795・MA_0826・MA_4115の3クエリを再実行し、以下を確認:

1. `06_Functional_Domain_Evidence`シートで、`radical_sam x electron_carrier`の一致候補が、修正前の32候補相当から、真のフラボドキシン3件(MA_0361, MA_0363, MA_2078)のみに絞り込まれていること。特に`MA_1812`(PF03358)・`MA_4170`(Pfamヒットなし)が一致しなくなっていることを確認。
2. MA_4115クエリでも同様に、`atp_dependent_activator x electron_carrier`の一致がMA_0361/MA_0363/MA_2078のみに絞られること(MA_0361/MA_0363自体は元々一致していたので、ここは非回帰確認)。
3. 全体の`final_score`ランキングがどう変化するか(電子伝達体候補が減ることで、他の候補が相対的にどう動くか)を確認。

## 受け入れ基準

- `python -m pytest -q`が全件パスすること。
- `electron_carrier`カテゴリの実データマッチ数が32件から3件に減ること。
- `cysteine_desulfurase`カテゴリの実データマッチ数が23件から6件に減ること。
- `dsre_sulfur_relay`カテゴリがMA_3928(DsrH)を含むようになり、かつMA_0809/MA_0398を含まなくなること。
- 既存の非回帰対象(MA_0361/MA_0363のatp_dependent_activator一致、MA_0331のiron_sulfur_cluster一致等)が壊れていないこと。

## 非破壊要件

- `category_rules`のカテゴリ名・ルール構成は変更しない(YAML内のアクセッションリストと`match_fields`のみの変更)。
- `domain_family_map_path`を設定していないユーザーへの影響はゼロ。
- `legacy_additive`スコアリングモデルの挙動には一切触れない。

## 対象外

- `cysteine_desulfurase`のPfamファミリー内部の粒度問題(PF00266がIscS/SufS型と一部アミノトランスフェラーゼを両方含む)の解消 — 配列レベルの追加判定が必要になるため、今回は見送り。
- 他生物種(archaea全般)での同カテゴリ定義の妥当性検証 — 今回はM. acetivorans実データのみでの監査。将来的にarchaea横断での汎用性を検証する場合は、他のarchaeaゲノムのUniProtバルクデータでも同じ監査(pfamのみ vs pfam+interpro+supfamの比較)を行うことを推奨する。

## 参考

- `claude/domain_family_map_v2_sulfur_relay_expansion.md`(v2設計・実データ検証ドキュメント。`match_fields`機構の元の実装経緯)
- `claude_code_instructions_domain_family_map_v2.md`(`match_fields`/`exclusive`のコード実装指示。今回のv3はコード変更なし、YAML更新のみ)
