
package us.kbase.kbhmmer;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: HMMER_MT_Bioelement_Params</p>
 * <pre>
 * HMMER MT_Bioelement Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_MT_Bioelement_N_ids",
    "input_MT_Bioelement_H_ids",
    "input_MT_Bioelement_O_ids",
    "input_MT_Bioelement_CFix_ids",
    "input_MT_Bioelement_C1_ids",
    "input_MT_Bioelement_CH4_ids",
    "input_MT_Bioelement_CO_ids",
    "input_MT_Bioelement_S_ids",
    "input_MT_Bioelement_CN_ids",
    "input_MT_Bioelement_CH4N2O_ids",
    "input_MT_Bioelement_Se_ids",
    "input_MT_Bioelement_Metal_ids",
    "input_MT_Bioelement_As_ids",
    "input_MT_Bioelement_Halo_ids",
    "input_many_refs",
    "output_filtered_name",
    "genome_disp_name_config",
    "count_category",
    "use_model_specific_thresholds",
    "show_target_block_headers",
    "coalesce_output",
    "save_ALL_featureSets",
    "save_ANY_featureSets",
    "e_value",
    "bitscore",
    "model_cov_perc",
    "maxaccepts",
    "heatmap",
    "low_val",
    "vertical",
    "show_blanks"
})
public class HMMERMTBioelementParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_MT_Bioelement_N_ids")
    private String inputMTBioelementNIds;
    @JsonProperty("input_MT_Bioelement_H_ids")
    private String inputMTBioelementHIds;
    @JsonProperty("input_MT_Bioelement_O_ids")
    private String inputMTBioelementOIds;
    @JsonProperty("input_MT_Bioelement_CFix_ids")
    private String inputMTBioelementCFixIds;
    @JsonProperty("input_MT_Bioelement_C1_ids")
    private String inputMTBioelementC1Ids;
    @JsonProperty("input_MT_Bioelement_CH4_ids")
    private String inputMTBioelementCH4Ids;
    @JsonProperty("input_MT_Bioelement_CO_ids")
    private String inputMTBioelementCOIds;
    @JsonProperty("input_MT_Bioelement_S_ids")
    private String inputMTBioelementSIds;
    @JsonProperty("input_MT_Bioelement_CN_ids")
    private String inputMTBioelementCNIds;
    @JsonProperty("input_MT_Bioelement_CH4N2O_ids")
    private String inputMTBioelementCH4N2OIds;
    @JsonProperty("input_MT_Bioelement_Se_ids")
    private String inputMTBioelementSeIds;
    @JsonProperty("input_MT_Bioelement_Metal_ids")
    private String inputMTBioelementMetalIds;
    @JsonProperty("input_MT_Bioelement_As_ids")
    private String inputMTBioelementAsIds;
    @JsonProperty("input_MT_Bioelement_Halo_ids")
    private String inputMTBioelementHaloIds;
    @JsonProperty("input_many_refs")
    private String inputManyRefs;
    @JsonProperty("output_filtered_name")
    private String outputFilteredName;
    @JsonProperty("genome_disp_name_config")
    private String genomeDispNameConfig;
    @JsonProperty("count_category")
    private String countCategory;
    @JsonProperty("use_model_specific_thresholds")
    private Long useModelSpecificThresholds;
    @JsonProperty("show_target_block_headers")
    private Long showTargetBlockHeaders;
    @JsonProperty("coalesce_output")
    private Long coalesceOutput;
    @JsonProperty("save_ALL_featureSets")
    private Long saveALLFeatureSets;
    @JsonProperty("save_ANY_featureSets")
    private Long saveANYFeatureSets;
    @JsonProperty("e_value")
    private Double eValue;
    @JsonProperty("bitscore")
    private Double bitscore;
    @JsonProperty("model_cov_perc")
    private Double modelCovPerc;
    @JsonProperty("maxaccepts")
    private Double maxaccepts;
    @JsonProperty("heatmap")
    private Long heatmap;
    @JsonProperty("low_val")
    private Long lowVal;
    @JsonProperty("vertical")
    private Long vertical;
    @JsonProperty("show_blanks")
    private Long showBlanks;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public HMMERMTBioelementParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_N_ids")
    public String getInputMTBioelementNIds() {
        return inputMTBioelementNIds;
    }

    @JsonProperty("input_MT_Bioelement_N_ids")
    public void setInputMTBioelementNIds(String inputMTBioelementNIds) {
        this.inputMTBioelementNIds = inputMTBioelementNIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementNIds(String inputMTBioelementNIds) {
        this.inputMTBioelementNIds = inputMTBioelementNIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_H_ids")
    public String getInputMTBioelementHIds() {
        return inputMTBioelementHIds;
    }

    @JsonProperty("input_MT_Bioelement_H_ids")
    public void setInputMTBioelementHIds(String inputMTBioelementHIds) {
        this.inputMTBioelementHIds = inputMTBioelementHIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementHIds(String inputMTBioelementHIds) {
        this.inputMTBioelementHIds = inputMTBioelementHIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_O_ids")
    public String getInputMTBioelementOIds() {
        return inputMTBioelementOIds;
    }

    @JsonProperty("input_MT_Bioelement_O_ids")
    public void setInputMTBioelementOIds(String inputMTBioelementOIds) {
        this.inputMTBioelementOIds = inputMTBioelementOIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementOIds(String inputMTBioelementOIds) {
        this.inputMTBioelementOIds = inputMTBioelementOIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_CFix_ids")
    public String getInputMTBioelementCFixIds() {
        return inputMTBioelementCFixIds;
    }

    @JsonProperty("input_MT_Bioelement_CFix_ids")
    public void setInputMTBioelementCFixIds(String inputMTBioelementCFixIds) {
        this.inputMTBioelementCFixIds = inputMTBioelementCFixIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementCFixIds(String inputMTBioelementCFixIds) {
        this.inputMTBioelementCFixIds = inputMTBioelementCFixIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_C1_ids")
    public String getInputMTBioelementC1Ids() {
        return inputMTBioelementC1Ids;
    }

    @JsonProperty("input_MT_Bioelement_C1_ids")
    public void setInputMTBioelementC1Ids(String inputMTBioelementC1Ids) {
        this.inputMTBioelementC1Ids = inputMTBioelementC1Ids;
    }

    public HMMERMTBioelementParams withInputMTBioelementC1Ids(String inputMTBioelementC1Ids) {
        this.inputMTBioelementC1Ids = inputMTBioelementC1Ids;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_CH4_ids")
    public String getInputMTBioelementCH4Ids() {
        return inputMTBioelementCH4Ids;
    }

    @JsonProperty("input_MT_Bioelement_CH4_ids")
    public void setInputMTBioelementCH4Ids(String inputMTBioelementCH4Ids) {
        this.inputMTBioelementCH4Ids = inputMTBioelementCH4Ids;
    }

    public HMMERMTBioelementParams withInputMTBioelementCH4Ids(String inputMTBioelementCH4Ids) {
        this.inputMTBioelementCH4Ids = inputMTBioelementCH4Ids;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_CO_ids")
    public String getInputMTBioelementCOIds() {
        return inputMTBioelementCOIds;
    }

    @JsonProperty("input_MT_Bioelement_CO_ids")
    public void setInputMTBioelementCOIds(String inputMTBioelementCOIds) {
        this.inputMTBioelementCOIds = inputMTBioelementCOIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementCOIds(String inputMTBioelementCOIds) {
        this.inputMTBioelementCOIds = inputMTBioelementCOIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_S_ids")
    public String getInputMTBioelementSIds() {
        return inputMTBioelementSIds;
    }

    @JsonProperty("input_MT_Bioelement_S_ids")
    public void setInputMTBioelementSIds(String inputMTBioelementSIds) {
        this.inputMTBioelementSIds = inputMTBioelementSIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementSIds(String inputMTBioelementSIds) {
        this.inputMTBioelementSIds = inputMTBioelementSIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_CN_ids")
    public String getInputMTBioelementCNIds() {
        return inputMTBioelementCNIds;
    }

    @JsonProperty("input_MT_Bioelement_CN_ids")
    public void setInputMTBioelementCNIds(String inputMTBioelementCNIds) {
        this.inputMTBioelementCNIds = inputMTBioelementCNIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementCNIds(String inputMTBioelementCNIds) {
        this.inputMTBioelementCNIds = inputMTBioelementCNIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_CH4N2O_ids")
    public String getInputMTBioelementCH4N2OIds() {
        return inputMTBioelementCH4N2OIds;
    }

    @JsonProperty("input_MT_Bioelement_CH4N2O_ids")
    public void setInputMTBioelementCH4N2OIds(String inputMTBioelementCH4N2OIds) {
        this.inputMTBioelementCH4N2OIds = inputMTBioelementCH4N2OIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementCH4N2OIds(String inputMTBioelementCH4N2OIds) {
        this.inputMTBioelementCH4N2OIds = inputMTBioelementCH4N2OIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_Se_ids")
    public String getInputMTBioelementSeIds() {
        return inputMTBioelementSeIds;
    }

    @JsonProperty("input_MT_Bioelement_Se_ids")
    public void setInputMTBioelementSeIds(String inputMTBioelementSeIds) {
        this.inputMTBioelementSeIds = inputMTBioelementSeIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementSeIds(String inputMTBioelementSeIds) {
        this.inputMTBioelementSeIds = inputMTBioelementSeIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_Metal_ids")
    public String getInputMTBioelementMetalIds() {
        return inputMTBioelementMetalIds;
    }

    @JsonProperty("input_MT_Bioelement_Metal_ids")
    public void setInputMTBioelementMetalIds(String inputMTBioelementMetalIds) {
        this.inputMTBioelementMetalIds = inputMTBioelementMetalIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementMetalIds(String inputMTBioelementMetalIds) {
        this.inputMTBioelementMetalIds = inputMTBioelementMetalIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_As_ids")
    public String getInputMTBioelementAsIds() {
        return inputMTBioelementAsIds;
    }

    @JsonProperty("input_MT_Bioelement_As_ids")
    public void setInputMTBioelementAsIds(String inputMTBioelementAsIds) {
        this.inputMTBioelementAsIds = inputMTBioelementAsIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementAsIds(String inputMTBioelementAsIds) {
        this.inputMTBioelementAsIds = inputMTBioelementAsIds;
        return this;
    }

    @JsonProperty("input_MT_Bioelement_Halo_ids")
    public String getInputMTBioelementHaloIds() {
        return inputMTBioelementHaloIds;
    }

    @JsonProperty("input_MT_Bioelement_Halo_ids")
    public void setInputMTBioelementHaloIds(String inputMTBioelementHaloIds) {
        this.inputMTBioelementHaloIds = inputMTBioelementHaloIds;
    }

    public HMMERMTBioelementParams withInputMTBioelementHaloIds(String inputMTBioelementHaloIds) {
        this.inputMTBioelementHaloIds = inputMTBioelementHaloIds;
        return this;
    }

    @JsonProperty("input_many_refs")
    public String getInputManyRefs() {
        return inputManyRefs;
    }

    @JsonProperty("input_many_refs")
    public void setInputManyRefs(String inputManyRefs) {
        this.inputManyRefs = inputManyRefs;
    }

    public HMMERMTBioelementParams withInputManyRefs(String inputManyRefs) {
        this.inputManyRefs = inputManyRefs;
        return this;
    }

    @JsonProperty("output_filtered_name")
    public String getOutputFilteredName() {
        return outputFilteredName;
    }

    @JsonProperty("output_filtered_name")
    public void setOutputFilteredName(String outputFilteredName) {
        this.outputFilteredName = outputFilteredName;
    }

    public HMMERMTBioelementParams withOutputFilteredName(String outputFilteredName) {
        this.outputFilteredName = outputFilteredName;
        return this;
    }

    @JsonProperty("genome_disp_name_config")
    public String getGenomeDispNameConfig() {
        return genomeDispNameConfig;
    }

    @JsonProperty("genome_disp_name_config")
    public void setGenomeDispNameConfig(String genomeDispNameConfig) {
        this.genomeDispNameConfig = genomeDispNameConfig;
    }

    public HMMERMTBioelementParams withGenomeDispNameConfig(String genomeDispNameConfig) {
        this.genomeDispNameConfig = genomeDispNameConfig;
        return this;
    }

    @JsonProperty("count_category")
    public String getCountCategory() {
        return countCategory;
    }

    @JsonProperty("count_category")
    public void setCountCategory(String countCategory) {
        this.countCategory = countCategory;
    }

    public HMMERMTBioelementParams withCountCategory(String countCategory) {
        this.countCategory = countCategory;
        return this;
    }

    @JsonProperty("use_model_specific_thresholds")
    public Long getUseModelSpecificThresholds() {
        return useModelSpecificThresholds;
    }

    @JsonProperty("use_model_specific_thresholds")
    public void setUseModelSpecificThresholds(Long useModelSpecificThresholds) {
        this.useModelSpecificThresholds = useModelSpecificThresholds;
    }

    public HMMERMTBioelementParams withUseModelSpecificThresholds(Long useModelSpecificThresholds) {
        this.useModelSpecificThresholds = useModelSpecificThresholds;
        return this;
    }

    @JsonProperty("show_target_block_headers")
    public Long getShowTargetBlockHeaders() {
        return showTargetBlockHeaders;
    }

    @JsonProperty("show_target_block_headers")
    public void setShowTargetBlockHeaders(Long showTargetBlockHeaders) {
        this.showTargetBlockHeaders = showTargetBlockHeaders;
    }

    public HMMERMTBioelementParams withShowTargetBlockHeaders(Long showTargetBlockHeaders) {
        this.showTargetBlockHeaders = showTargetBlockHeaders;
        return this;
    }

    @JsonProperty("coalesce_output")
    public Long getCoalesceOutput() {
        return coalesceOutput;
    }

    @JsonProperty("coalesce_output")
    public void setCoalesceOutput(Long coalesceOutput) {
        this.coalesceOutput = coalesceOutput;
    }

    public HMMERMTBioelementParams withCoalesceOutput(Long coalesceOutput) {
        this.coalesceOutput = coalesceOutput;
        return this;
    }

    @JsonProperty("save_ALL_featureSets")
    public Long getSaveALLFeatureSets() {
        return saveALLFeatureSets;
    }

    @JsonProperty("save_ALL_featureSets")
    public void setSaveALLFeatureSets(Long saveALLFeatureSets) {
        this.saveALLFeatureSets = saveALLFeatureSets;
    }

    public HMMERMTBioelementParams withSaveALLFeatureSets(Long saveALLFeatureSets) {
        this.saveALLFeatureSets = saveALLFeatureSets;
        return this;
    }

    @JsonProperty("save_ANY_featureSets")
    public Long getSaveANYFeatureSets() {
        return saveANYFeatureSets;
    }

    @JsonProperty("save_ANY_featureSets")
    public void setSaveANYFeatureSets(Long saveANYFeatureSets) {
        this.saveANYFeatureSets = saveANYFeatureSets;
    }

    public HMMERMTBioelementParams withSaveANYFeatureSets(Long saveANYFeatureSets) {
        this.saveANYFeatureSets = saveANYFeatureSets;
        return this;
    }

    @JsonProperty("e_value")
    public Double getEValue() {
        return eValue;
    }

    @JsonProperty("e_value")
    public void setEValue(Double eValue) {
        this.eValue = eValue;
    }

    public HMMERMTBioelementParams withEValue(Double eValue) {
        this.eValue = eValue;
        return this;
    }

    @JsonProperty("bitscore")
    public Double getBitscore() {
        return bitscore;
    }

    @JsonProperty("bitscore")
    public void setBitscore(Double bitscore) {
        this.bitscore = bitscore;
    }

    public HMMERMTBioelementParams withBitscore(Double bitscore) {
        this.bitscore = bitscore;
        return this;
    }

    @JsonProperty("model_cov_perc")
    public Double getModelCovPerc() {
        return modelCovPerc;
    }

    @JsonProperty("model_cov_perc")
    public void setModelCovPerc(Double modelCovPerc) {
        this.modelCovPerc = modelCovPerc;
    }

    public HMMERMTBioelementParams withModelCovPerc(Double modelCovPerc) {
        this.modelCovPerc = modelCovPerc;
        return this;
    }

    @JsonProperty("maxaccepts")
    public Double getMaxaccepts() {
        return maxaccepts;
    }

    @JsonProperty("maxaccepts")
    public void setMaxaccepts(Double maxaccepts) {
        this.maxaccepts = maxaccepts;
    }

    public HMMERMTBioelementParams withMaxaccepts(Double maxaccepts) {
        this.maxaccepts = maxaccepts;
        return this;
    }

    @JsonProperty("heatmap")
    public Long getHeatmap() {
        return heatmap;
    }

    @JsonProperty("heatmap")
    public void setHeatmap(Long heatmap) {
        this.heatmap = heatmap;
    }

    public HMMERMTBioelementParams withHeatmap(Long heatmap) {
        this.heatmap = heatmap;
        return this;
    }

    @JsonProperty("low_val")
    public Long getLowVal() {
        return lowVal;
    }

    @JsonProperty("low_val")
    public void setLowVal(Long lowVal) {
        this.lowVal = lowVal;
    }

    public HMMERMTBioelementParams withLowVal(Long lowVal) {
        this.lowVal = lowVal;
        return this;
    }

    @JsonProperty("vertical")
    public Long getVertical() {
        return vertical;
    }

    @JsonProperty("vertical")
    public void setVertical(Long vertical) {
        this.vertical = vertical;
    }

    public HMMERMTBioelementParams withVertical(Long vertical) {
        this.vertical = vertical;
        return this;
    }

    @JsonProperty("show_blanks")
    public Long getShowBlanks() {
        return showBlanks;
    }

    @JsonProperty("show_blanks")
    public void setShowBlanks(Long showBlanks) {
        this.showBlanks = showBlanks;
    }

    public HMMERMTBioelementParams withShowBlanks(Long showBlanks) {
        this.showBlanks = showBlanks;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("HMMERMTBioelementParams"+" [workspaceName=")+ workspaceName)+", inputMTBioelementNIds=")+ inputMTBioelementNIds)+", inputMTBioelementHIds=")+ inputMTBioelementHIds)+", inputMTBioelementOIds=")+ inputMTBioelementOIds)+", inputMTBioelementCFixIds=")+ inputMTBioelementCFixIds)+", inputMTBioelementC1Ids=")+ inputMTBioelementC1Ids)+", inputMTBioelementCH4Ids=")+ inputMTBioelementCH4Ids)+", inputMTBioelementCOIds=")+ inputMTBioelementCOIds)+", inputMTBioelementSIds=")+ inputMTBioelementSIds)+", inputMTBioelementCNIds=")+ inputMTBioelementCNIds)+", inputMTBioelementCH4N2OIds=")+ inputMTBioelementCH4N2OIds)+", inputMTBioelementSeIds=")+ inputMTBioelementSeIds)+", inputMTBioelementMetalIds=")+ inputMTBioelementMetalIds)+", inputMTBioelementAsIds=")+ inputMTBioelementAsIds)+", inputMTBioelementHaloIds=")+ inputMTBioelementHaloIds)+", inputManyRefs=")+ inputManyRefs)+", outputFilteredName=")+ outputFilteredName)+", genomeDispNameConfig=")+ genomeDispNameConfig)+", countCategory=")+ countCategory)+", useModelSpecificThresholds=")+ useModelSpecificThresholds)+", showTargetBlockHeaders=")+ showTargetBlockHeaders)+", coalesceOutput=")+ coalesceOutput)+", saveALLFeatureSets=")+ saveALLFeatureSets)+", saveANYFeatureSets=")+ saveANYFeatureSets)+", eValue=")+ eValue)+", bitscore=")+ bitscore)+", modelCovPerc=")+ modelCovPerc)+", maxaccepts=")+ maxaccepts)+", heatmap=")+ heatmap)+", lowVal=")+ lowVal)+", vertical=")+ vertical)+", showBlanks=")+ showBlanks)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
