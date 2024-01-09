package ru.biosoft.physicell.core;

@FunctionalInterface
public interface PhaseEntry
{
    public void execute(Cell pCell, Phenotype phenotype, double dt);

    public static class standard_Ki67_positive_phase_entry_function implements PhaseEntry
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // the cell wants to double its volume 
            phenotype.volume.target_solid_nuclear *= 2.0;
            phenotype.volume.target_solid_cytoplasmic *= 2.0;
        }
    }

    public static class standard_Ki67_negative_phase_entry_function implements PhaseEntry
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            return;
        }
    }

    public static class StandardLivePhaseEntry implements PhaseEntry
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // the cell wants to double its volume 
            phenotype.volume.target_solid_nuclear *= 2.0;
            phenotype.volume.target_solid_cytoplasmic *= 2.0;
        }
    }

    public static class SPhaseEntry implements PhaseEntry
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // the cell wants to double its volume 
            phenotype.volume.target_solid_nuclear *= 2.0;
            phenotype.volume.target_solid_cytoplasmic *= 2.0;
        }
    }

    public static class StandardCyclingEntry implements PhaseEntry
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // the cell wants to double its volume 
            phenotype.volume.target_solid_nuclear *= 2.0;
            phenotype.volume.target_solid_cytoplasmic *= 2.0;
        }
    }

    public static class standard_apoptosis_entry_function implements PhaseEntry
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // the volume model wants to shrink the cell
            Volume v = phenotype.volume;
            v.target_fluid_fraction = 0.0;
            v.target_solid_cytoplasmic = 0.0;
            v.target_solid_nuclear = 0.0;

            v.target_cytoplasmic_to_nuclear_ratio = 0.0;

            // change the rate parameters 
            DeathParameters p = phenotype.death.current_parameters();
            v.cytoplasmic_biomass_change_rate = p.cytoplasmic_biomass_change_rate;
            v.nuclear_biomass_change_rate = p.nuclear_biomass_change_rate;
            v.fluid_change_rate = p.unlysed_fluid_change_rate;

            v.calcification_rate = p.calcification_rate;

            v.relative_rupture_volume = p.relative_rupture_volume;
            v.rupture_volume = v.total * v.relative_rupture_volume;
        }
    }

    public static class standard_necrosis_entry_function implements PhaseEntry
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // the volume model wants to degrade the solids, but swell by osmosis 
            Volume v = phenotype.volume;
            v.target_fluid_fraction = 1.0;
            v.target_solid_cytoplasmic = 0.0;
            v.target_solid_nuclear = 0.0;

            v.target_cytoplasmic_to_nuclear_ratio = 0.0;

            // change the rate parameters   
            DeathParameters p = phenotype.death.current_parameters();
            v.cytoplasmic_biomass_change_rate = p.cytoplasmic_biomass_change_rate;
            v.nuclear_biomass_change_rate = p.nuclear_biomass_change_rate;
            v.fluid_change_rate = p.unlysed_fluid_change_rate;

            v.calcification_rate = p.calcification_rate;

            // set the bursting volume 
            v.relative_rupture_volume = p.relative_rupture_volume;
            v.rupture_volume = v.total * v.relative_rupture_volume;
        }
    }

    public static class standard_lysis_entry_function implements PhaseEntry
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // the volume model wants to shrink the cell
            Volume v = phenotype.volume;
            v.target_fluid_fraction = 0.0;
            v.target_solid_cytoplasmic = 0.0;
            v.target_solid_nuclear = 0.0;

            // change the rate parameters 
            DeathParameters p = phenotype.death.current_parameters();
            v.cytoplasmic_biomass_change_rate = p.cytoplasmic_biomass_change_rate;
            v.nuclear_biomass_change_rate = p.nuclear_biomass_change_rate;
            v.fluid_change_rate = p.lysed_fluid_change_rate;

            v.calcification_rate = p.calcification_rate;

            // set the bursting volume 
            v.relative_rupture_volume = 9e99;
            v.rupture_volume = v.total * v.relative_rupture_volume;
        }
    }
}