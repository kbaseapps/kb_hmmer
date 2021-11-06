
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
 * <p>Original spec-file type: HMMER_EnvBioelement_Params</p>
 * <pre>
 * HMMER EnvBioelement Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_EnvBioelement_N_ids",
    "input_EnvBioelement_H_ids",
    "input_EnvBioelement_O_ids",
    "input_EnvBioelement_CFix_ids",
    "input_EnvBioelement_C1_ids",
    "input_EnvBioelement_CH4_ids",
    "input_EnvBioelement_CO_ids",
    "input_EnvBioelement_S_ids",
    "input_EnvBioelement_CN_ids",
    "input_EnvBioelement_CH4N2O_ids",
    "input_EnvBioelement_Se_ids",
    "input_EnvBioelement_Metal_ids",
    "input_EnvBioelement_As_ids",
    "input_EnvBioelement_Halo_ids",
    "input_many_refs",
    "output_filtered_name",
    "genome_disp_name_config",
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
public class HMMEREnvBioelementParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_EnvBioelement_N_ids")
    private String inputEnvBioelementNIds;
    @JsonProperty("input_EnvBioelement_H_ids")
    private String inputEnvBioelementHIds;
    @JsonProperty("input_EnvBioelement_O_ids")
    private String inputEnvBioelementOIds;
    @JsonProperty("input_EnvBioelement_CFix_ids")
    private String inputEnvBioelementCFixIds;
    @JsonProperty("input_EnvBioelement_C1_ids")
    private String inputEnvBioelementC1Ids;
    @JsonProperty("input_EnvBioelement_CH4_ids")
    private String inputEnvBioelementCH4Ids;
    @JsonProperty("input_EnvBioelement_CO_ids")
    private String inputEnvBioelementCOIds;
    @JsonProperty("input_EnvBioelement_S_ids")
    private String inputEnvBioelementSIds;
    @JsonProperty("input_EnvBioelement_CN_ids")
    private String inputEnvBioelementCNIds;
    @JsonProperty("input_EnvBioelement_CH4N2O_ids")
    private String inputEnvBioelementCH4N2OIds;
    @JsonProperty("input_EnvBioelement_Se_ids")
    private String inputEnvBioelementSeIds;
    @JsonProperty("input_EnvBioelement_Metal_ids")
    private String inputEnvBioelementMetalIds;
    @JsonProperty("input_EnvBioelement_As_ids")
    private String inputEnvBioelementAsIds;
    @JsonProperty("input_EnvBioelement_Halo_ids")
    private String inputEnvBioelementHaloIds;
    @JsonProperty("input_many_refs")
    private String inputManyRefs;
    @JsonProperty("output_filtered_name")
    private String outputFilteredName;
    @JsonProperty("genome_disp_name_config")
    private String genomeDispNameConfig;
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

    public HMMEREnvBioelementParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_EnvBioelement_N_ids")
    public String getInputEnvBioelementNIds() {
        return inputEnvBioelementNIds;
    }

    @JsonProperty("input_EnvBioelement_N_ids")
    public void setInputEnvBioelementNIds(String inputEnvBioelementNIds) {
        this.inputEnvBioelementNIds = inputEnvBioelementNIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementNIds(String inputEnvBioelementNIds) {
        this.inputEnvBioelementNIds = inputEnvBioelementNIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_H_ids")
    public String getInputEnvBioelementHIds() {
        return inputEnvBioelementHIds;
    }

    @JsonProperty("input_EnvBioelement_H_ids")
    public void setInputEnvBioelementHIds(String inputEnvBioelementHIds) {
        this.inputEnvBioelementHIds = inputEnvBioelementHIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementHIds(String inputEnvBioelementHIds) {
        this.inputEnvBioelementHIds = inputEnvBioelementHIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_O_ids")
    public String getInputEnvBioelementOIds() {
        return inputEnvBioelementOIds;
    }

    @JsonProperty("input_EnvBioelement_O_ids")
    public void setInputEnvBioelementOIds(String inputEnvBioelementOIds) {
        this.inputEnvBioelementOIds = inputEnvBioelementOIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementOIds(String inputEnvBioelementOIds) {
        this.inputEnvBioelementOIds = inputEnvBioelementOIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_CFix_ids")
    public String getInputEnvBioelementCFixIds() {
        return inputEnvBioelementCFixIds;
    }

    @JsonProperty("input_EnvBioelement_CFix_ids")
    public void setInputEnvBioelementCFixIds(String inputEnvBioelementCFixIds) {
        this.inputEnvBioelementCFixIds = inputEnvBioelementCFixIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementCFixIds(String inputEnvBioelementCFixIds) {
        this.inputEnvBioelementCFixIds = inputEnvBioelementCFixIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_C1_ids")
    public String getInputEnvBioelementC1Ids() {
        return inputEnvBioelementC1Ids;
    }

    @JsonProperty("input_EnvBioelement_C1_ids")
    public void setInputEnvBioelementC1Ids(String inputEnvBioelementC1Ids) {
        this.inputEnvBioelementC1Ids = inputEnvBioelementC1Ids;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementC1Ids(String inputEnvBioelementC1Ids) {
        this.inputEnvBioelementC1Ids = inputEnvBioelementC1Ids;
        return this;
    }

    @JsonProperty("input_EnvBioelement_CH4_ids")
    public String getInputEnvBioelementCH4Ids() {
        return inputEnvBioelementCH4Ids;
    }

    @JsonProperty("input_EnvBioelement_CH4_ids")
    public void setInputEnvBioelementCH4Ids(String inputEnvBioelementCH4Ids) {
        this.inputEnvBioelementCH4Ids = inputEnvBioelementCH4Ids;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementCH4Ids(String inputEnvBioelementCH4Ids) {
        this.inputEnvBioelementCH4Ids = inputEnvBioelementCH4Ids;
        return this;
    }

    @JsonProperty("input_EnvBioelement_CO_ids")
    public String getInputEnvBioelementCOIds() {
        return inputEnvBioelementCOIds;
    }

    @JsonProperty("input_EnvBioelement_CO_ids")
    public void setInputEnvBioelementCOIds(String inputEnvBioelementCOIds) {
        this.inputEnvBioelementCOIds = inputEnvBioelementCOIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementCOIds(String inputEnvBioelementCOIds) {
        this.inputEnvBioelementCOIds = inputEnvBioelementCOIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_S_ids")
    public String getInputEnvBioelementSIds() {
        return inputEnvBioelementSIds;
    }

    @JsonProperty("input_EnvBioelement_S_ids")
    public void setInputEnvBioelementSIds(String inputEnvBioelementSIds) {
        this.inputEnvBioelementSIds = inputEnvBioelementSIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementSIds(String inputEnvBioelementSIds) {
        this.inputEnvBioelementSIds = inputEnvBioelementSIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_CN_ids")
    public String getInputEnvBioelementCNIds() {
        return inputEnvBioelementCNIds;
    }

    @JsonProperty("input_EnvBioelement_CN_ids")
    public void setInputEnvBioelementCNIds(String inputEnvBioelementCNIds) {
        this.inputEnvBioelementCNIds = inputEnvBioelementCNIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementCNIds(String inputEnvBioelementCNIds) {
        this.inputEnvBioelementCNIds = inputEnvBioelementCNIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_CH4N2O_ids")
    public String getInputEnvBioelementCH4N2OIds() {
        return inputEnvBioelementCH4N2OIds;
    }

    @JsonProperty("input_EnvBioelement_CH4N2O_ids")
    public void setInputEnvBioelementCH4N2OIds(String inputEnvBioelementCH4N2OIds) {
        this.inputEnvBioelementCH4N2OIds = inputEnvBioelementCH4N2OIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementCH4N2OIds(String inputEnvBioelementCH4N2OIds) {
        this.inputEnvBioelementCH4N2OIds = inputEnvBioelementCH4N2OIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_Se_ids")
    public String getInputEnvBioelementSeIds() {
        return inputEnvBioelementSeIds;
    }

    @JsonProperty("input_EnvBioelement_Se_ids")
    public void setInputEnvBioelementSeIds(String inputEnvBioelementSeIds) {
        this.inputEnvBioelementSeIds = inputEnvBioelementSeIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementSeIds(String inputEnvBioelementSeIds) {
        this.inputEnvBioelementSeIds = inputEnvBioelementSeIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_Metal_ids")
    public String getInputEnvBioelementMetalIds() {
        return inputEnvBioelementMetalIds;
    }

    @JsonProperty("input_EnvBioelement_Metal_ids")
    public void setInputEnvBioelementMetalIds(String inputEnvBioelementMetalIds) {
        this.inputEnvBioelementMetalIds = inputEnvBioelementMetalIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementMetalIds(String inputEnvBioelementMetalIds) {
        this.inputEnvBioelementMetalIds = inputEnvBioelementMetalIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_As_ids")
    public String getInputEnvBioelementAsIds() {
        return inputEnvBioelementAsIds;
    }

    @JsonProperty("input_EnvBioelement_As_ids")
    public void setInputEnvBioelementAsIds(String inputEnvBioelementAsIds) {
        this.inputEnvBioelementAsIds = inputEnvBioelementAsIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementAsIds(String inputEnvBioelementAsIds) {
        this.inputEnvBioelementAsIds = inputEnvBioelementAsIds;
        return this;
    }

    @JsonProperty("input_EnvBioelement_Halo_ids")
    public String getInputEnvBioelementHaloIds() {
        return inputEnvBioelementHaloIds;
    }

    @JsonProperty("input_EnvBioelement_Halo_ids")
    public void setInputEnvBioelementHaloIds(String inputEnvBioelementHaloIds) {
        this.inputEnvBioelementHaloIds = inputEnvBioelementHaloIds;
    }

    public HMMEREnvBioelementParams withInputEnvBioelementHaloIds(String inputEnvBioelementHaloIds) {
        this.inputEnvBioelementHaloIds = inputEnvBioelementHaloIds;
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

    public HMMEREnvBioelementParams withInputManyRefs(String inputManyRefs) {
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

    public HMMEREnvBioelementParams withOutputFilteredName(String outputFilteredName) {
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

    public HMMEREnvBioelementParams withGenomeDispNameConfig(String genomeDispNameConfig) {
        this.genomeDispNameConfig = genomeDispNameConfig;
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

    public HMMEREnvBioelementParams withShowTargetBlockHeaders(Long showTargetBlockHeaders) {
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

    public HMMEREnvBioelementParams withCoalesceOutput(Long coalesceOutput) {
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

    public HMMEREnvBioelementParams withSaveALLFeatureSets(Long saveALLFeatureSets) {
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

    public HMMEREnvBioelementParams withSaveANYFeatureSets(Long saveANYFeatureSets) {
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

    public HMMEREnvBioelementParams withEValue(Double eValue) {
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

    public HMMEREnvBioelementParams withBitscore(Double bitscore) {
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

    public HMMEREnvBioelementParams withModelCovPerc(Double modelCovPerc) {
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

    public HMMEREnvBioelementParams withMaxaccepts(Double maxaccepts) {
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

    public HMMEREnvBioelementParams withHeatmap(Long heatmap) {
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

    public HMMEREnvBioelementParams withLowVal(Long lowVal) {
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

    public HMMEREnvBioelementParams withVertical(Long vertical) {
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

    public HMMEREnvBioelementParams withShowBlanks(Long showBlanks) {
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
        return ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("HMMEREnvBioelementParams"+" [workspaceName=")+ workspaceName)+", inputEnvBioelementNIds=")+ inputEnvBioelementNIds)+", inputEnvBioelementHIds=")+ inputEnvBioelementHIds)+", inputEnvBioelementOIds=")+ inputEnvBioelementOIds)+", inputEnvBioelementCFixIds=")+ inputEnvBioelementCFixIds)+", inputEnvBioelementC1Ids=")+ inputEnvBioelementC1Ids)+", inputEnvBioelementCH4Ids=")+ inputEnvBioelementCH4Ids)+", inputEnvBioelementCOIds=")+ inputEnvBioelementCOIds)+", inputEnvBioelementSIds=")+ inputEnvBioelementSIds)+", inputEnvBioelementCNIds=")+ inputEnvBioelementCNIds)+", inputEnvBioelementCH4N2OIds=")+ inputEnvBioelementCH4N2OIds)+", inputEnvBioelementSeIds=")+ inputEnvBioelementSeIds)+", inputEnvBioelementMetalIds=")+ inputEnvBioelementMetalIds)+", inputEnvBioelementAsIds=")+ inputEnvBioelementAsIds)+", inputEnvBioelementHaloIds=")+ inputEnvBioelementHaloIds)+", inputManyRefs=")+ inputManyRefs)+", outputFilteredName=")+ outputFilteredName)+", genomeDispNameConfig=")+ genomeDispNameConfig)+", showTargetBlockHeaders=")+ showTargetBlockHeaders)+", coalesceOutput=")+ coalesceOutput)+", saveALLFeatureSets=")+ saveALLFeatureSets)+", saveANYFeatureSets=")+ saveANYFeatureSets)+", eValue=")+ eValue)+", bitscore=")+ bitscore)+", modelCovPerc=")+ modelCovPerc)+", maxaccepts=")+ maxaccepts)+", heatmap=")+ heatmap)+", lowVal=")+ lowVal)+", vertical=")+ vertical)+", showBlanks=")+ showBlanks)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
