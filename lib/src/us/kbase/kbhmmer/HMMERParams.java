
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
 * <p>Original spec-file type: HMMER_Params</p>
 * <pre>
 * HMMER Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_many_name",
    "input_msa_name",
    "output_filtered_name"
})
public class HMMERParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_many_name")
    private String inputManyName;
    @JsonProperty("input_msa_name")
    private String inputMsaName;
    @JsonProperty("output_filtered_name")
    private String outputFilteredName;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public HMMERParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("input_many_name")
    public String getInputManyName() {
        return inputManyName;
    }

    @JsonProperty("input_many_name")
    public void setInputManyName(String inputManyName) {
        this.inputManyName = inputManyName;
    }

    public HMMERParams withInputManyName(String inputManyName) {
        this.inputManyName = inputManyName;
        return this;
    }

    @JsonProperty("input_msa_name")
    public String getInputMsaName() {
        return inputMsaName;
    }

    @JsonProperty("input_msa_name")
    public void setInputMsaName(String inputMsaName) {
        this.inputMsaName = inputMsaName;
    }

    public HMMERParams withInputMsaName(String inputMsaName) {
        this.inputMsaName = inputMsaName;
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

    public HMMERParams withOutputFilteredName(String outputFilteredName) {
        this.outputFilteredName = outputFilteredName;
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
        return ((((((((((("HMMERParams"+" [workspaceName=")+ workspaceName)+", inputManyName=")+ inputManyName)+", inputMsaName=")+ inputMsaName)+", outputFilteredName=")+ outputFilteredName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
