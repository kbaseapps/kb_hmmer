
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
 * <p>Original spec-file type: HMMER_PhyloMarkers_Params</p>
 * <pre>
 * HMMER PhyloMarkers Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_PhyloMarkers_Univ_ids",
    "input_PhyloMarkers_B_ribo_pol_ids",
    "input_PhyloMarkers_B_other_ids",
    "input_PhyloMarkers_A_ribo_pol_ids",
    "input_PhyloMarkers_A_other_ids",
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
public class HMMERPhyloMarkersParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_PhyloMarkers_Univ_ids")
    private String inputPhyloMarkersUnivIds;
    @JsonProperty("input_PhyloMarkers_B_ribo_pol_ids")
    private String inputPhyloMarkersBRiboPolIds;
    @JsonProperty("input_PhyloMarkers_B_other_ids")
    private String inputPhyloMarkersBOtherIds;
    @JsonProperty("input_PhyloMarkers_A_ribo_pol_ids")
    private String inputPhyloMarkersARiboPolIds;
    @JsonProperty("input_PhyloMarkers_A_other_ids")
    private String inputPhyloMarkersAOtherIds;
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

    public HMMERPhyloMarkersParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_PhyloMarkers_Univ_ids")
    public String getInputPhyloMarkersUnivIds() {
        return inputPhyloMarkersUnivIds;
    }

    @JsonProperty("input_PhyloMarkers_Univ_ids")
    public void setInputPhyloMarkersUnivIds(String inputPhyloMarkersUnivIds) {
        this.inputPhyloMarkersUnivIds = inputPhyloMarkersUnivIds;
    }

    public HMMERPhyloMarkersParams withInputPhyloMarkersUnivIds(String inputPhyloMarkersUnivIds) {
        this.inputPhyloMarkersUnivIds = inputPhyloMarkersUnivIds;
        return this;
    }

    @JsonProperty("input_PhyloMarkers_B_ribo_pol_ids")
    public String getInputPhyloMarkersBRiboPolIds() {
        return inputPhyloMarkersBRiboPolIds;
    }

    @JsonProperty("input_PhyloMarkers_B_ribo_pol_ids")
    public void setInputPhyloMarkersBRiboPolIds(String inputPhyloMarkersBRiboPolIds) {
        this.inputPhyloMarkersBRiboPolIds = inputPhyloMarkersBRiboPolIds;
    }

    public HMMERPhyloMarkersParams withInputPhyloMarkersBRiboPolIds(String inputPhyloMarkersBRiboPolIds) {
        this.inputPhyloMarkersBRiboPolIds = inputPhyloMarkersBRiboPolIds;
        return this;
    }

    @JsonProperty("input_PhyloMarkers_B_other_ids")
    public String getInputPhyloMarkersBOtherIds() {
        return inputPhyloMarkersBOtherIds;
    }

    @JsonProperty("input_PhyloMarkers_B_other_ids")
    public void setInputPhyloMarkersBOtherIds(String inputPhyloMarkersBOtherIds) {
        this.inputPhyloMarkersBOtherIds = inputPhyloMarkersBOtherIds;
    }

    public HMMERPhyloMarkersParams withInputPhyloMarkersBOtherIds(String inputPhyloMarkersBOtherIds) {
        this.inputPhyloMarkersBOtherIds = inputPhyloMarkersBOtherIds;
        return this;
    }

    @JsonProperty("input_PhyloMarkers_A_ribo_pol_ids")
    public String getInputPhyloMarkersARiboPolIds() {
        return inputPhyloMarkersARiboPolIds;
    }

    @JsonProperty("input_PhyloMarkers_A_ribo_pol_ids")
    public void setInputPhyloMarkersARiboPolIds(String inputPhyloMarkersARiboPolIds) {
        this.inputPhyloMarkersARiboPolIds = inputPhyloMarkersARiboPolIds;
    }

    public HMMERPhyloMarkersParams withInputPhyloMarkersARiboPolIds(String inputPhyloMarkersARiboPolIds) {
        this.inputPhyloMarkersARiboPolIds = inputPhyloMarkersARiboPolIds;
        return this;
    }

    @JsonProperty("input_PhyloMarkers_A_other_ids")
    public String getInputPhyloMarkersAOtherIds() {
        return inputPhyloMarkersAOtherIds;
    }

    @JsonProperty("input_PhyloMarkers_A_other_ids")
    public void setInputPhyloMarkersAOtherIds(String inputPhyloMarkersAOtherIds) {
        this.inputPhyloMarkersAOtherIds = inputPhyloMarkersAOtherIds;
    }

    public HMMERPhyloMarkersParams withInputPhyloMarkersAOtherIds(String inputPhyloMarkersAOtherIds) {
        this.inputPhyloMarkersAOtherIds = inputPhyloMarkersAOtherIds;
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

    public HMMERPhyloMarkersParams withInputManyRefs(String inputManyRefs) {
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

    public HMMERPhyloMarkersParams withOutputFilteredName(String outputFilteredName) {
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

    public HMMERPhyloMarkersParams withGenomeDispNameConfig(String genomeDispNameConfig) {
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

    public HMMERPhyloMarkersParams withShowTargetBlockHeaders(Long showTargetBlockHeaders) {
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

    public HMMERPhyloMarkersParams withCoalesceOutput(Long coalesceOutput) {
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

    public HMMERPhyloMarkersParams withSaveALLFeatureSets(Long saveALLFeatureSets) {
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

    public HMMERPhyloMarkersParams withSaveANYFeatureSets(Long saveANYFeatureSets) {
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

    public HMMERPhyloMarkersParams withEValue(Double eValue) {
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

    public HMMERPhyloMarkersParams withBitscore(Double bitscore) {
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

    public HMMERPhyloMarkersParams withModelCovPerc(Double modelCovPerc) {
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

    public HMMERPhyloMarkersParams withMaxaccepts(Double maxaccepts) {
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

    public HMMERPhyloMarkersParams withHeatmap(Long heatmap) {
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

    public HMMERPhyloMarkersParams withLowVal(Long lowVal) {
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

    public HMMERPhyloMarkersParams withVertical(Long vertical) {
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

    public HMMERPhyloMarkersParams withShowBlanks(Long showBlanks) {
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
        return ((((((((((((((((((((((((((((((((((((((((((((("HMMERPhyloMarkersParams"+" [workspaceName=")+ workspaceName)+", inputPhyloMarkersUnivIds=")+ inputPhyloMarkersUnivIds)+", inputPhyloMarkersBRiboPolIds=")+ inputPhyloMarkersBRiboPolIds)+", inputPhyloMarkersBOtherIds=")+ inputPhyloMarkersBOtherIds)+", inputPhyloMarkersARiboPolIds=")+ inputPhyloMarkersARiboPolIds)+", inputPhyloMarkersAOtherIds=")+ inputPhyloMarkersAOtherIds)+", inputManyRefs=")+ inputManyRefs)+", outputFilteredName=")+ outputFilteredName)+", genomeDispNameConfig=")+ genomeDispNameConfig)+", showTargetBlockHeaders=")+ showTargetBlockHeaders)+", coalesceOutput=")+ coalesceOutput)+", saveALLFeatureSets=")+ saveALLFeatureSets)+", saveANYFeatureSets=")+ saveANYFeatureSets)+", eValue=")+ eValue)+", bitscore=")+ bitscore)+", modelCovPerc=")+ modelCovPerc)+", maxaccepts=")+ maxaccepts)+", heatmap=")+ heatmap)+", lowVal=")+ lowVal)+", vertical=")+ vertical)+", showBlanks=")+ showBlanks)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
