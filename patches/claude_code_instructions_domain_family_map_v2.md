# Claude Code 実装指示書: domain_family_map v2 (硫黄リレー・tRNAチオール化ドメイン拡張)

## 背景・目的

`domain_complementarity` v3 (`analysis/domain_family_map.py`, コミット `81b4c44`) は MA_4115 実データで検証済み (`claude/domain_complementarity_v3_design.md` 参照)。ただし現在 `config/domain_family_map.v1.yaml` に実データが入っているカテゴリは `electron_carrier` と `atp_dependent_activator` の2つのみで、`sulfur_carrier` と `radical_sam` は空 (アクセッション未設定) のプレースホルダのままです。

今回、ユーザーの研究テーマ (*Methanosarcina acetivorans* の cnm⁵U 生合成経路・MA_4115 下流探索) に直結する硫黄リレー/tRNAチオール化系のドメインファミリーについて、UniProt バルクエクスポート実データ (`uniprotkb_methanosarcina_acetivorans_2026_09_06.json`) と、他生物種の構造生物学文献 (ThiS-ThiF: PDB 1ZUD、IscS-IscU: PLOS Biology 2010、TusA-TusBCD-TusE: Structure 2006 ほか) を調査し、新規5カテゴリ・6ルールの追加と、既存2カテゴリの修正を設計しました。詳細な調査経緯・文献根拠・実データ表は `claude/domain_family_map_v2_sulfur_relay_expansion.md` を参照してください (このセッションでユーザー承認済み)。

このセッション (Cowork) は設計・調査のみを行っており、実際のコード変更は行っていません。以下の指示に従って `~/projects/ProteinHunter` に実装してください。

## 変更概要

| 変更 | 内容 |
|---|---|
| コード変更 | `analysis/domain_family_map.py`: `DomainFamilyCategory` に `exclusive` / `match_fields` の2フィールドを追加し、`categories_for()` に排他的一致・フィールド限定一致のロジックを追加 |
| 設定変更 | `config/domain_family_map.v1.yaml`: カテゴリを4→9個、`category_rules` を3→9件に拡張 |
| テスト追加 | `tests/test_domain_family_map.py` に `exclusive`/`match_fields` の単体テストを追加、`tests/test_interaction_scoring.py` に統合テストを1〜2件追加 (推奨) |

**この変更は `config/domain_family_map.v1.yaml` を明示的に設定しているユーザーにのみ影響します。** `domain_family_map_path` が未設定のユーザーには一切影響しません (`_domain_complementarity_status_and_value` の既存の `None` チェックによりレイヤーごとスキップされる、既存の非破壊設計はそのまま)。

## 対象ファイル

- `analysis/domain_family_map.py` (コード変更)
- `config/domain_family_map.v1.yaml` (設定の全面差し替え)
- `tests/test_domain_family_map.py` (テスト追加)
- `tests/test_interaction_scoring.py` (テスト追加、推奨)

## 実装詳細

### 1. `analysis/domain_family_map.py`: `exclusive` / `match_fields` の追加

現状の `categories_for()` は「候補の Pfam/InterPro/SUPFAM のいずれかがカテゴリのアクセッションと**交差**すれば一致」という部分一致ロジックのみです。今回、2つの新しい一致モードが必要になります。

#### 1-1. なぜ必要か

- **`exclusive: true` (新設、`iron_sulfur_cluster` カテゴリで使用)**: PF00037 (Fer4) は単独フェレドキシンだけでなく、ヘテロジスルフィド還元酵素・ABC輸送体・グルタミン酸合成酵素など大きな多ドメイン複合体サブユニットにも広く出現します (実データで計48件ヒット)。単独フェレドキシン (MA_1485 など、Pfam ヒットが `[PF00037]` のみ) だけを `iron_sulfur_cluster` として扱いたいので、「候補の Pfam ヒット**全体**がカテゴリの Pfam アクセッションの部分集合である場合にのみ一致」という排他的一致モードが必要です。
- **`match_fields` (新設、`atp_dependent_activator` カテゴリで使用)**: 現在の `atp_dependent_activator` は `interpro: [IPR014729]` / `supfam: [SSF52402]` という HUP スーパーファミリーの構造フォールド全般を指すアクセッションを含んでおり、これは MA_4115 だけでなく無関係な ATP/AMP 利用酵素群 (例えば MA_1466 / thiI) にも共通するフォールドです。実際 MA_1466 は `IPR014729` と `SSF52402` を MA_4115 と共有しているため、現状のロジックのままだと `atp_dependent_activator` が意図せず広く一致してしまいます。`match_fields: [pfam]` を指定できるようにし、このカテゴリだけ Pfam (`PF24167`) のみで判定するようにします。interpro/supfam のアクセッション自体は参考情報として YAML に残しますが、一致判定には使いません。

#### 1-2. `DomainFamilyCategory` データクラスの変更

```python
@dataclass(frozen=True, slots=True)
class DomainFamilyCategory:
    """One domain-family category's Pfam/InterPro/SUPFAM accession set."""

    pfam: frozenset[str] = frozenset()
    interpro: frozenset[str] = frozenset()
    supfam: frozenset[str] = frozenset()
    exclusive: bool = False
    match_fields: frozenset[str] = frozenset({"pfam", "interpro", "supfam"})
```

`exclusive`/`match_fields` を追加してもデフォルト値により既存の呼び出し・既存テストは変更不要です (後方互換)。

#### 1-3. `categories_for()` の変更

```python
def categories_for(self, info: "UniProtDomainInfo | None") -> set[str]:
    """Return every category ``info`` belongs to (empty set if ``info`` is None)."""
    if info is None:
        return set()

    matched: set[str] = set()
    for category_id, category in self.categories.items():
        if category.exclusive:
            # Exclusive mode: match only when the candidate's ENTIRE pfam hit
            # set is a subset of this category's pfam accessions (i.e. the
            # candidate carries no other, unrelated Pfam domains). An empty
            # candidate pfam set must NOT match -- the empty set is
            # mathematically a subset of any set, but a protein with zero
            # Pfam hits carries no positive evidence of belonging here.
            if info.pfam and info.pfam <= category.pfam:
                matched.add(category_id)
            continue

        hit = (
            ("pfam" in category.match_fields and bool(info.pfam & category.pfam))
            or ("interpro" in category.match_fields and bool(info.interpro & category.interpro))
            or ("supfam" in category.match_fields and bool(info.supfam & category.supfam))
        )
        if hit:
            matched.add(category_id)
    return matched
```

`find_category_match()` は変更不要です (カテゴリ集合同士の比較のみで、一致モードの詳細を知る必要がないため)。

#### 1-4. `load_domain_family_map()` の変更 (YAMLパース)

`categories` のループ内、`DomainFamilyCategory(...)` を構築する直前に以下を追加してください:

```python
        raw_exclusive = raw_category.get("exclusive", False)
        if not isinstance(raw_exclusive, bool):
            raise ConfigError(
                f"'categories.{category_id}.exclusive' in {resolved_path} must be true or false."
            )

        raw_match_fields = raw_category.get("match_fields")
        if raw_match_fields is None:
            match_fields = frozenset({"pfam", "interpro", "supfam"})
        else:
            match_fields = frozenset(
                _string_list(raw_match_fields, f"categories.{category_id}.match_fields", resolved_path)
            )
            unknown_fields = match_fields - {"pfam", "interpro", "supfam"}
            if unknown_fields:
                raise ConfigError(
                    f"'categories.{category_id}.match_fields' in {resolved_path} "
                    f"contains unknown field(s): {sorted(unknown_fields)} "
                    "(allowed: pfam, interpro, supfam)."
                )
            if not match_fields:
                raise ConfigError(
                    f"'categories.{category_id}.match_fields' in {resolved_path} must not be empty."
                )
```

そして既存の `categories[str(category_id)] = DomainFamilyCategory(...)` の呼び出しに `exclusive=raw_exclusive, match_fields=match_fields` を追加してください。

`exclusive: true` かつ `match_fields` が明示指定されているケースは今回のカテゴリ案には存在しませんが (`iron_sulfur_cluster` は `exclusive: true` のみ、`match_fields` はデフォルトのまま使わない — `exclusive` ブランチでは `match_fields` を参照しないため無害です)、将来のため両立可能な設計にしてあります。

### 2. `config/domain_family_map.v1.yaml` の全面差し替え

現在のファイルには `electron_carrier` (実データ入り)・`sulfur_carrier` (空)・`atp_dependent_activator` (実データ入り)・`radical_sam` (空、note のみ) の4カテゴリと、3件の `category_rules` があります。**「`radical_sam` カテゴリが未定義」という問題ではなく、カテゴリ自体は既にありますが Pfam/InterPro/SUPFAM が空リストのため `categories_for()` が絶対にヒットしない (空集合との積は必ず空) という状態です。** 以下の内容で `config/domain_family_map.v1.yaml` 全体を置き換えてください (既存の `electron_carrier`/`atp_dependent_activator` の実データ・note は保持し、`sulfur_carrier`/`radical_sam` を実データで埋め、新規5カテゴリ・6ルールを追加する形です):

```yaml
version: v1
categories:
  electron_carrier:
    pfam: [PF12724]
    interpro: [IPR026816, IPR029039]
    supfam: [SSF52218]
    note: "flavodoxin/flavoprotein様の電子伝達体。MA_0361/MA_0363で確認済み(2026-09-09)。"

  atp_dependent_activator:
    pfam: [PF24167]
    interpro: [IPR055834, IPR014729]   # 参考情報のみ。match_fieldsによりpfamのみで判定。
    supfam: [SSF52402]                  # 同上。
    match_fields: [pfam]
    note: >
      HUP/PP-loop ATPaseドメイン。MA_4115で確認済み(2026-09-09)。
      IPR014729/SSF52402はHUPスーパーファミリー全般(ThiI等の無関係な酵素も含む広い構造フォールド)を
      指すため判定には使わず、match_fields: [pfam] でPF24167のみによる判定に限定する(2026-09-09改訂)。

  sulfur_carrier:
    pfam: [PF02597, PF21965]   # ThiS/MoaD(PF02597), SAMP2(PF21965)
    interpro: [IPR003749, IPR052045]   # ThiS/MoaD-like, Sulfur_Carrier/Prot_Modifier
    supfam: [SSF54285]   # MoaD/ThiS
    note: "ubiquitin-fold型の硫黄キャリア/タンパク質修飾因子。実データ: MA_4086(moaD), MA_1713(MoaD様), MA_3300(SAMP2)。(2026-09-09追加)"

  cysteine_desulfurase:
    pfam: [PF00266]   # Aminotran_5
    interpro: [IPR010240, IPR010970, IPR016454]   # Cys_deSase_IscS, Cys_dSase_SufS, Cysteine_dSase
    supfam: [SSF53383]   # PLP-dependent transferases
    note: "システイン脱硫酵素(IscS/SufS型)。実データ: MA_3264(iscS), MA_2718(iscS), MA_0236, MA_0808。(2026-09-09新規)"

  rhodanese_sulfurtransferase:
    pfam: [PF00581]
    interpro: [IPR001763]
    supfam: [SSF52821]
    note: "ロダネーゼ様硫黄転移ドメイン。実データ: MA_3178, MA_2500, MA_0746。(2026-09-09新規)"

  thif_moeb_activator:
    pfam: [PF00899]   # ThiF
    interpro: [IPR045886, IPR035985]   # ThiF/MoeB/HesA, Ubiquitin-activating_enz
    supfam: [SSF69572]
    note: "E1様活性化酵素(ThiS/MoaD/SAMPをアデニル化する側)。実データ: MA_0255。(2026-09-09新規)"

  tusa_sulfur_relay:
    pfam: [PF01206]   # TusA
    interpro: [IPR001455]   # TusA-like
    note: >
      TusA型硫黄リレータンパク質(直接パースルフィド中継、非アデニル化)。
      実データ: MA_0810(現状「UPF0033 domain-containing protein」、機能未確定)。(2026-09-09新規)

  dsre_sulfur_relay:
    pfam: [PF02635]   # DsrE
    interpro: [IPR027396, IPR003787]   # DsrEFH-like, Sulphur_relay_DsrE/F-like
    note: >
      TusBCD/DsrEFH型硫黄リレータンパク質(TusAから硫黄を受け取る側)。
      実データ: MA_3926, MA_3927(隣接遺伝子、オペロンの可能性), MA_3072。いずれも現状「Uncharacterized protein」、機能未確定。(2026-09-09新規)

  radical_sam:
    pfam: [PF04055]
    note: >
      radical SAM酵素。実データ40件ヒット、うちtRNA修飾関連: MA_0040(taw1), MA_1153(t(6)A37 methylthiotransferase),
      MA_0826(hat)。既存ルールradical_sam×electron_carrierが参照していたが、これまで本カテゴリが空リストのため
      一致しなかった不具合を修正(2026-09-09)。

  iron_sulfur_cluster:
    pfam: [PF00037]   # Fer4
    exclusive: true
    note: >
      単独フェレドキシン限定(排他的一致)。実データ: MA_1485(Pfamヒットが[PF00037]のみ)は該当。
      PF00037系Fer4ドメインはヘテロジスルフィド還元酵素・ABC輸送体・グルタミン酸合成酵素など大きな
      多ドメイン複合体サブユニットにも広く出現する(計48件ヒット)ため、exclusive: trueにより、
      Pfamヒット全体が[PF00037]の部分集合である候補のみに限定する。(2026-09-09新規)

category_rules:
  - left: atp_dependent_activator
    right: electron_carrier
    note: "ATP依存アデニル化酵素(HUP/PP-loopドメイン)と電子伝達体/フラボプロテイン様タンパク質"
  - left: atp_dependent_activator
    right: sulfur_carrier
    note: "ATP依存アデニル化酵素(HUP/PP-loopドメイン)と硫黄キャリア"
  - left: radical_sam
    right: electron_carrier
    note: "radical SAM酵素と電子伝達体(2026-09-09: radical_samカテゴリに実データを設定し機能するように修正)"
  - left: radical_sam
    right: iron_sulfur_cluster
    note: "radical SAM酵素(自身も[4Fe-4S]クラスタを持つ)とFe-S電子伝達体(2026-09-09新規)"
  - left: cysteine_desulfurase
    right: sulfur_carrier
    note: "システイン脱硫酵素(硫黄供与)と硫黄キャリアタンパク質(ThiS/MoaD/SAMP型)。IscS-IscU/SufS-SufE型相互作用の一般化(PLOS Biology 2010)。(2026-09-09新規)"
  - left: cysteine_desulfurase
    right: rhodanese_sulfurtransferase
    note: "システイン脱硫酵素とロダネーゼ様硫黄転移タンパク質(持続的硫黄リレー)。(2026-09-09新規)"
  - left: sulfur_carrier
    right: thif_moeb_activator
    note: "硫黄キャリア(ThiS/MoaD/SAMP型)とそれを活性化するE1様酵素(ThiF/MoeB型)。ThiS-ThiF複合体構造(PDB 1ZUD)、Urm1-Uba4系(EMBO J 2020)の一般化。(2026-09-09新規)"
  - left: tusa_sulfur_relay
    right: dsre_sulfur_relay
    note: >
      TusA型硫黄リレータンパク質とTusBCD/DsrEFH型受容タンパク質。大腸菌TusA-TusBCD-TusE系の複合体構造
      (Structure誌2006)の一般化。M. acetivoransにもMA_0810(TusA型)とMA_3926/3927/3072(DsrE型)の
      ホモログが存在するが、いずれも機能未確定。(2026-09-09新規)
  - left: cysteine_desulfurase
    right: tusa_sulfur_relay
    note: "システイン脱硫酵素とTusA型硫黄リレータンパク質(IscSからTusAへの硫黄受け渡し)。(2026-09-09新規)"
```

## テスト計画

### `tests/test_domain_family_map.py` (必須)

以下を新規テストとして追加してください (既存の `SAMPLE_YAML` に依存するテストはそのままで問題ありません。`exclusive`/`match_fields` 用に別の小さな YAML 断片を各テスト内で用意してください):

1. **`test_exclusive_category_matches_only_when_pfam_is_a_subset`**: `exclusive: true` のカテゴリ (pfam: `[PF00037]`) に対し、`UniProtDomainInfo(pfam=frozenset({"PF00037"}))` (単独) は一致し、`UniProtDomainInfo(pfam=frozenset({"PF00037", "PF99999"}))` (他のPfamも持つ) は一致しないことを確認。
2. **`test_exclusive_category_does_not_match_empty_pfam`**: `UniProtDomainInfo(pfam=frozenset())` (Pfamヒットなし) が `exclusive: true` カテゴリに一致しないことを確認 (空集合はどんな集合の部分集合でもあるため、この保護が無いと誤って一致してしまう回帰を防ぐためのテスト)。
3. **`test_match_fields_restricts_matching_to_listed_fields`**: `match_fields: [pfam]` を指定したカテゴリに対し、interpro のみが一致する `UniProtDomainInfo` (pfamは不一致) が一致しないこと、pfam が一致する場合は一致することを確認。
4. **`test_shipped_v1_domain_family_map_loads` の更新**: 新しいカテゴリ数 (9) ・`category_rules` 件数 (9) を確認するアサーションを追加し、`loaded.categories["iron_sulfur_cluster"].exclusive is True`、`loaded.categories["atp_dependent_activator"].match_fields == frozenset({"pfam"})`、`loaded.categories["radical_sam"].pfam == frozenset({"PF04055"})` を確認。

### `tests/test_interaction_scoring.py` (推奨)

既存の `test_v2_mode_domain_family_map_scores_pfam_match_without_legacy_keyword_overlap` と同じパターンで、以下を追加することを推奨します (必須ではありません):

5. **`test_v2_mode_domain_family_map_exclusive_category_excludes_multidomain_candidate`**: `iron_sulfur_cluster` (exclusive) 相当のテスト用ドメインマップを用意し、Fer4 ドメイン単体を持つ候補ではカテゴリマッチが成立し、Fer4 + 別の無関係ドメインを持つ候補では成立しない (= `domain_complementarity` がヒットしない) ことをエンドツーエンドで確認。

## 動作確認手順

```bash
cd ~/projects/ProteinHunter
python -m pytest tests/test_domain_family_map.py -q
python -m pytest -q   # 全件(テスト追加分だけ既存の523から増加するはず)
```

その後、MA_4115 クエリで実際にパイプラインを再実行し、以下を確認してください:

1. `06_Functional_Domain_Evidence` シートで、新たに `explanation` に `domain family match (UniProt Pfam/InterPro): cysteine_desulfurase x sulfur_carrier` 等、今回追加したカテゴリペアの一致が現れる候補が (もしあれば) 出てくること。
2. MA_1466 (thiI) が `atp_dependent_activator` として MA_4115 と一致しなくなっている (match_fields: [pfam] 適用により、以前 interpro/supfam 経由で誤って一致していた可能性がある場合の回帰確認) こと — ただし現時点でこの誤マッチが実際に発生していたかどうかは未確認のため、「一致しないこと」自体を新たな成功基準にするのではなく、「意図通りpfam限定で判定されているか `explanation` を確認する」という位置づけで構いません。
3. 既存の MA_0361/MA_0363 (`atp_dependent_activator x electron_carrier`) の一致が壊れていないこと (非回帰確認)。

## 受け入れ基準

- `python -m pytest -q` が全件パスすること。
- `config/domain_family_map.v1.yaml` が `load_domain_family_map()` でエラーなくロードできること (`test_shipped_v1_domain_family_map_loads` が該当)。
- `exclusive: true` カテゴリが、単独ドメイン候補にのみ一致し、多ドメイン候補には一致しないこと (単体テストで確認)。
- `match_fields: [pfam]` を指定した `atp_dependent_activator` が、Pfam不一致・InterPro/SUPFAMのみ一致のケースで一致しなくなること (単体テストで確認)。
- 既存の `electron_carrier`/`atp_dependent_activator` の実データ一致 (MA_0361/MA_0363, MA_4115) が壊れていないこと。

## 非破壊要件

- `domain_family_map_path` を設定していないユーザーへの影響はゼロ (既存の `None` チェックによりレイヤーごとスキップされる設計は変更しない)。
- `legacy_additive` スコアリングモデルの挙動には一切触れない。
- `exclusive`/`match_fields` はどちらもデフォルト値 (`False` / 全フィールド) を持つため、既存の `categories:` エントリ (`exclusive`/`match_fields` を指定していないもの) の挙動は完全に従来通り。

## 対象外

- MA_1466 (thiI)・MA_3300 (SAMP2)・MA_0810/MA_3926系列 (TusA/DsrE型) を MA_4115 の直接相互作用候補として個別に文献調査すること — ユーザー承認により今回は見送り。
- IscU/NifU型 Fe-Sクラスタ scaffold タンパク質、SufE の Pfam アクセッション確定 — 今回のUniProtバルクデータ調査で確証が得られなかったため、カテゴリに含めていない (`claude/domain_family_map_v2_sulfur_relay_expansion.md` の「未確定・今回見送ったもの」参照)。

## 参考

- `claude/domain_family_map_v2_sulfur_relay_expansion.md` (このセッションでの調査・設計ドキュメント本体。文献根拠・実データ表を含む)
- `claude/domain_complementarity_v3_design.md` (v3の元設計・実データ検証結果)
