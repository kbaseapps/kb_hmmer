
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
 * <p>Original spec-file type: HMMER_Local_MSA_Group_Params</p>
 * <pre>
 * HMMER Local MSA Group Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_msa_refs",
    "input_many_ref",
    "output_filtered_name",
    "genome_disp_name_config",
    "coalesce_output",
    "e_value",
    "bitscore",
    "model_cov_perc",
    "maxaccepts",
    "heatmap",
    "low_val",
    "vertical",
    "show_blanks"
})
public class HMMERLocalMSAGroupParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_msa_refs")
    private String inputMsaRefs;
    @JsonProperty("input_many_ref")
    private String inputManyRef;
    @JsonProperty("output_filtered_name")
    private String outputFilteredName;
    @JsonProperty("genome_disp_name_config")
    private String genomeDispNameConfig;
    @JsonProperty("coalesce_output")
    private Long coalesceOutput;
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

    public HMMERLocalMSAGroupParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_msa_refs")
    public String getInputMsaRefs() {
        return inputMsaRefs;
    }

    @JsonProperty("input_msa_refs")
    public void setInputMsaRefs(String inputMsaRefs) {
        this.inputMsaRefs = inputMsaRefs;
    }

    public HMMERLocalMSAGroupParams withInputMsaRefs(String inputMsaRefs) {
        this.inputMsaRefs = inputMsaRefs;
        return this;
    }

    @JsonProperty("input_many_ref")
    public String getInputManyRef() {
        return inputManyRef;
    }

    @JsonProperty("input_many_ref")
    public void setInputManyRef(String inputManyRef) {
        this.inputManyRef = inputManyRef;
    }

    public HMMERLocalMSAGroupParams withInputManyRef(String inputManyRef) {
        this.inputManyRef = inputManyRef;
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

    public HMMERLocalMSAGroupParams withOutputFilteredName(String outputFilteredName) {
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

    public HMMERLocalMSAGroupParams withGenomeDispNameConfig(String genomeDispNameConfig) {
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

    public HMMERLocalMSAGroupParams withCoalesceOutput(Long coalesceOutput) {
        this.coalesceOutput = coalesceOutput;
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

    public HMMERLocalMSAGroupParams withEValue(Double eValue) {
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

    public HMMERLocalMSAGroupParams withBitscore(Double bitscore) {
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

    public HMMERLocalMSAGroupParams withModelCovPerc(Double modelCovPerc) {
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

    public HMMERLocalMSAGroupParams withMaxaccepts(Double maxaccepts) {
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

    public HMMERLocalMSAGroupParams withHeatmap(Long heatmap) {
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

    public HMMERLocalMSAGroupParams withLowVal(Long lowVal) {
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

    public HMMERLocalMSAGroupParams withVertical(Long vertical) {
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

    public HMMERLocalMSAGroupParams withShowBlanks(Long showBlanks) {
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
        return ((((((((((((((((((((((((((((((("HMMERLocalMSAGroupParams"+" [workspaceName=")+ workspaceName)+", inputMsaRefs=")+ inputMsaRefs)+", inputManyRef=")+ inputManyRef)+", outputFilteredName=")+ outputFilteredName)+", genomeDispNameConfig=")+ genomeDispNameConfig)+", coalesceOutput=")+ coalesceOutput)+", eValue=")+ eValue)+", bitscore=")+ bitscore)+", modelCovPerc=")+ modelCovPerc)+", maxaccepts=")+ maxaccepts)+", heatmap=")+ heatmap)+", lowVal=")+ lowVal)+", vertical=")+ vertical)+", showBlanks=")+ showBlanks)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
