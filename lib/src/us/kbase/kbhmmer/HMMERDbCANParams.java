
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
 * <p>Original spec-file type: HMMER_dbCAN_Params</p>
 * <pre>
 * HMMER dbCAN Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_dbCAN_AA_ids",
    "input_dbCAN_CBM_ids",
    "input_dbCAN_CE_ids",
    "input_dbCAN_GH_ids",
    "input_dbCAN_GT_ids",
    "input_dbCAN_PL_ids",
    "input_dbCAN_cellulosome_ids",
    "input_many_refs",
    "output_filtered_name",
    "genome_disp_name_config",
    "coalesce_output",
    "save_ALL_featureSets",
    "e_value",
    "bitscore",
    "model_cov_perc",
    "maxaccepts",
    "heatmap",
    "low_val",
    "vertical",
    "show_blanks"
})
public class HMMERDbCANParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_dbCAN_AA_ids")
    private String inputDbCANAAIds;
    @JsonProperty("input_dbCAN_CBM_ids")
    private String inputDbCANCBMIds;
    @JsonProperty("input_dbCAN_CE_ids")
    private String inputDbCANCEIds;
    @JsonProperty("input_dbCAN_GH_ids")
    private String inputDbCANGHIds;
    @JsonProperty("input_dbCAN_GT_ids")
    private String inputDbCANGTIds;
    @JsonProperty("input_dbCAN_PL_ids")
    private String inputDbCANPLIds;
    @JsonProperty("input_dbCAN_cellulosome_ids")
    private String inputDbCANCellulosomeIds;
    @JsonProperty("input_many_refs")
    private String inputManyRefs;
    @JsonProperty("output_filtered_name")
    private String outputFilteredName;
    @JsonProperty("genome_disp_name_config")
    private String genomeDispNameConfig;
    @JsonProperty("coalesce_output")
    private Long coalesceOutput;
    @JsonProperty("save_ALL_featureSets")
    private Long saveALLFeatureSets;
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

    public HMMERDbCANParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_dbCAN_AA_ids")
    public String getInputDbCANAAIds() {
        return inputDbCANAAIds;
    }

    @JsonProperty("input_dbCAN_AA_ids")
    public void setInputDbCANAAIds(String inputDbCANAAIds) {
        this.inputDbCANAAIds = inputDbCANAAIds;
    }

    public HMMERDbCANParams withInputDbCANAAIds(String inputDbCANAAIds) {
        this.inputDbCANAAIds = inputDbCANAAIds;
        return this;
    }

    @JsonProperty("input_dbCAN_CBM_ids")
    public String getInputDbCANCBMIds() {
        return inputDbCANCBMIds;
    }

    @JsonProperty("input_dbCAN_CBM_ids")
    public void setInputDbCANCBMIds(String inputDbCANCBMIds) {
        this.inputDbCANCBMIds = inputDbCANCBMIds;
    }

    public HMMERDbCANParams withInputDbCANCBMIds(String inputDbCANCBMIds) {
        this.inputDbCANCBMIds = inputDbCANCBMIds;
        return this;
    }

    @JsonProperty("input_dbCAN_CE_ids")
    public String getInputDbCANCEIds() {
        return inputDbCANCEIds;
    }

    @JsonProperty("input_dbCAN_CE_ids")
    public void setInputDbCANCEIds(String inputDbCANCEIds) {
        this.inputDbCANCEIds = inputDbCANCEIds;
    }

    public HMMERDbCANParams withInputDbCANCEIds(String inputDbCANCEIds) {
        this.inputDbCANCEIds = inputDbCANCEIds;
        return this;
    }

    @JsonProperty("input_dbCAN_GH_ids")
    public String getInputDbCANGHIds() {
        return inputDbCANGHIds;
    }

    @JsonProperty("input_dbCAN_GH_ids")
    public void setInputDbCANGHIds(String inputDbCANGHIds) {
        this.inputDbCANGHIds = inputDbCANGHIds;
    }

    public HMMERDbCANParams withInputDbCANGHIds(String inputDbCANGHIds) {
        this.inputDbCANGHIds = inputDbCANGHIds;
        return this;
    }

    @JsonProperty("input_dbCAN_GT_ids")
    public String getInputDbCANGTIds() {
        return inputDbCANGTIds;
    }

    @JsonProperty("input_dbCAN_GT_ids")
    public void setInputDbCANGTIds(String inputDbCANGTIds) {
        this.inputDbCANGTIds = inputDbCANGTIds;
    }

    public HMMERDbCANParams withInputDbCANGTIds(String inputDbCANGTIds) {
        this.inputDbCANGTIds = inputDbCANGTIds;
        return this;
    }

    @JsonProperty("input_dbCAN_PL_ids")
    public String getInputDbCANPLIds() {
        return inputDbCANPLIds;
    }

    @JsonProperty("input_dbCAN_PL_ids")
    public void setInputDbCANPLIds(String inputDbCANPLIds) {
        this.inputDbCANPLIds = inputDbCANPLIds;
    }

    public HMMERDbCANParams withInputDbCANPLIds(String inputDbCANPLIds) {
        this.inputDbCANPLIds = inputDbCANPLIds;
        return this;
    }

    @JsonProperty("input_dbCAN_cellulosome_ids")
    public String getInputDbCANCellulosomeIds() {
        return inputDbCANCellulosomeIds;
    }

    @JsonProperty("input_dbCAN_cellulosome_ids")
    public void setInputDbCANCellulosomeIds(String inputDbCANCellulosomeIds) {
        this.inputDbCANCellulosomeIds = inputDbCANCellulosomeIds;
    }

    public HMMERDbCANParams withInputDbCANCellulosomeIds(String inputDbCANCellulosomeIds) {
        this.inputDbCANCellulosomeIds = inputDbCANCellulosomeIds;
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

    public HMMERDbCANParams withInputManyRefs(String inputManyRefs) {
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

    public HMMERDbCANParams withOutputFilteredName(String outputFilteredName) {
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

    public HMMERDbCANParams withGenomeDispNameConfig(String genomeDispNameConfig) {
        this.genomeDispNameConfig = genomeDispNameConfig;
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

    public HMMERDbCANParams withCoalesceOutput(Long coalesceOutput) {
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

    public HMMERDbCANParams withSaveALLFeatureSets(Long saveALLFeatureSets) {
        this.saveALLFeatureSets = saveALLFeatureSets;
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

    public HMMERDbCANParams withEValue(Double eValue) {
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

    public HMMERDbCANParams withBitscore(Double bitscore) {
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

    public HMMERDbCANParams withModelCovPerc(Double modelCovPerc) {
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

    public HMMERDbCANParams withMaxaccepts(Double maxaccepts) {
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

    public HMMERDbCANParams withHeatmap(Long heatmap) {
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

    public HMMERDbCANParams withLowVal(Long lowVal) {
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

    public HMMERDbCANParams withVertical(Long vertical) {
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

    public HMMERDbCANParams withShowBlanks(Long showBlanks) {
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
        return ((((((((((((((((((((((((((((((((((((((((((((("HMMERDbCANParams"+" [workspaceName=")+ workspaceName)+", inputDbCANAAIds=")+ inputDbCANAAIds)+", inputDbCANCBMIds=")+ inputDbCANCBMIds)+", inputDbCANCEIds=")+ inputDbCANCEIds)+", inputDbCANGHIds=")+ inputDbCANGHIds)+", inputDbCANGTIds=")+ inputDbCANGTIds)+", inputDbCANPLIds=")+ inputDbCANPLIds)+", inputDbCANCellulosomeIds=")+ inputDbCANCellulosomeIds)+", inputManyRefs=")+ inputManyRefs)+", outputFilteredName=")+ outputFilteredName)+", genomeDispNameConfig=")+ genomeDispNameConfig)+", coalesceOutput=")+ coalesceOutput)+", saveALLFeatureSets=")+ saveALLFeatureSets)+", eValue=")+ eValue)+", bitscore=")+ bitscore)+", modelCovPerc=")+ modelCovPerc)+", maxaccepts=")+ maxaccepts)+", heatmap=")+ heatmap)+", lowVal=")+ lowVal)+", vertical=")+ vertical)+", showBlanks=")+ showBlanks)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
