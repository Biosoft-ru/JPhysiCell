package ru.biosoft.physicell.core;

public class PhysiCellSettings
{
    double max_time = 60 * 24 * 45;

    // units
    String time_units = "min";
    String space_units = "micron";

    // parallel options 
    int omp_num_threads = 2;

    // save options
    String folder = ".";

    double full_save_interval = 60;
    boolean enable_full_saves = true;
    boolean enable_legacy_saves = false;

    boolean disable_automated_spring_adhesions = false;

    double SVG_save_interval = 60;
    boolean enable_SVG_saves = true;

    boolean enable_substrate_plot = false;
    String substrate_to_monitor = "oxygen";
    boolean limits_substrate_plot = false;
    double min_concentration = -1.0;
    double max_concentration = -1.0;

    double intracellular_save_interval = 60;
    boolean enable_intracellular_saves = false;

    // cell rules option
    boolean rules_enabled = false;
    String rules_protocol = "Cell Behavior Hypothesis Grammar (CBHG)";
    String rules_protocol_version = "1.0";

    public PhysiCellSettings()
    {
        // units 
        time_units = "min";
        space_units = "micron";

        // save options
        folder = ".";
        max_time = 60 * 24 * 45;

        full_save_interval = 60;
        enable_full_saves = true;
        enable_legacy_saves = false;

        SVG_save_interval = 60;
        enable_SVG_saves = true;
        enable_substrate_plot = false;
        substrate_to_monitor = "oxygen";
        limits_substrate_plot = false;
        min_concentration = -1.0;
        max_concentration = -1.0;

        intracellular_save_interval = 60;
        enable_intracellular_saves = false;

        // parallel options 

        omp_num_threads = 4;

        rules_enabled = false;
        rules_protocol = "Cell Behavior Hypothesis Grammar (CBHG)";
        rules_protocol_version = "1.0";

        return;
    }

}
