# Audit: how COSMIC evidence can reach the main classification (OS4/OM4 leak)

Scope: every path by which `OS4`, `OM4`, `COSMIC_TOTAL_OCCURRENCES`, `COSMIC_MATCH`
or `cosmic_count` can influence `final_class`, `somatic_oncogenicity_class`,
`CONFIDENCE_SCORE_BASE`, `CONFIDENCE_CATEGORY_BASE`, `filter_status`,
`filter_reasons`. Verified by reading the code.

## Findings

| Regra | Arquivo/função | Entrada usada | Considera tipo tumoral? | Resultado afetado | Comportamento (após correção) |
| --- | --- | --- | --- | --- | --- |
| OS4 | `guideline_classification.R::add_somatic_oncogenicity_classification` (was L50, uses `cosmic_count`) | **antes:** `COSMIC_TOTAL_OCCURRENCES` (global) | **antes: NÃO** | `somatic_oncogenicity_score/class` → `final_class`/`filter_status` | **CORRIGIDO:** agora usa `OS4_CONTEXTUAL` (só ocorrências no tumor da amostra, `>= os4_min`) |
| OM4 | `guideline_classification.R::add_somatic_oncogenicity_classification` (was L51) | **antes:** `COSMIC_TOTAL_OCCURRENCES` (global) | **antes: NÃO** | `somatic_oncogenicity_score/class` → `final_class`/`filter_status` | **CORRIGIDO:** agora usa `OM4_CONTEXTUAL` (matched, `>= om4_min` e `< os4_min`) |
| `somatic_oncogenicity_class ∈ {Oncogenic, Likely Oncogenic}` | `scoring.R::classify_variants` L175-176 (`guideline_oncogenic`) | pontos de oncogenicidade (agora incluindo OS4/OM4 **contextuais**) | Sim (via OS4/OM4 contextual) | `final_class = high_confidence_somatic` | OK: só recebe OS4/OM4 quando o contexto tumoral permite |
| `oncogenic_category = "oncogenic_supportive"` (de `COSMIC_MATCH`) | `scoring.R::classify_oncogenic_evidence` L296 | `COSMIC_MATCH` (booleano) | Não | **NENHUM** — `classify_variants` usa apenas `oncogenic_exact` (L176/L195); `supportive` não é lido | Sem vazamento (confirmado) |
| `validation_support_score` (= `COSMIC_CONTEXT_SUPPORT_SCORE`) | `scoring.R::classify_variants` L127-128 | contexto tumoral (banded) | Sim | `somatic_score_validated` → **apenas** confiança EXPERIMENTAL | Não afeta `final_class` nem `CONFIDENCE_SCORE_BASE` |
| `cosmic_driver_score` (de `COSMIC_MATCH`/`cosmic_count` global) | `driver_classification.R::cosmic_driver_score` L243-268 | `COSMIC_TOTAL_OCCURRENCES` (global) + `COSMIC_MATCH` | Não | `driver_score`/`driver_class` (anotação separada) | **Não** afeta `final_class`/`filter_status`/`CONFIDENCE_SCORE_BASE`; é uma anotação de driver à parte (documentado) |
| `cosmic_count` como coluna | `cosmic_match.R::annotate_cosmic` | `COSMIC_TOTAL_OCCURRENCES` | — | usado só por driver (acima) | mantido como evidência global informativa |

## Conclusão

O único caminho pelo qual a recorrência **global** do COSMIC alterava a
classificação principal era OS4/OM4 em `guideline_classification.R`. Após a
correção, OS4/OM4 são atribuídos exclusivamente a partir de
`OS4_CONTEXTUAL`/`OM4_CONTEXTUAL`, que exigem ocorrências no tipo tumoral da
amostra (`COSMIC_MATCHING_TUMOR_OCCURRENCES`) e `COSMIC_CONTEXT_EVALUABLE == TRUE`.
Recorrência exclusiva em outros tumores (`other_tumor_only`,
`pan_cancer_without_sample_tumor`) e contextos não avaliáveis (`tumor_type_unknown`,
`cosmic_tumor_type_missing`) são marcados `COSMIC_GLOBAL_EVIDENCE_ONLY = TRUE` /
`COSMIC_CONTEXT_INTERPRETATION = not_evaluable` e **não** geram OS4/OM4.

`CONFIDENCE_SCORE_BASE`/`CONFIDENCE_CATEGORY_BASE` e `filter_status` permanecem
independentes de qualquer evidência COSMIC contextual (regressão testada).
