package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.DeathParameters;
import ru.biosoft.physicell.core.PhaseEntry;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.Volume;

public class standard_apoptosis_entry_function implements PhaseEntry
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
        DeathParameters p = phenotype.death.currentParameters();
        v.cytoplasmic_biomass_change_rate = p.cytoplasmic_biomass_change_rate;
        v.nuclear_biomass_change_rate = p.nuclear_biomass_change_rate;
        v.fluid_change_rate = p.unlysed_fluid_change_rate;

        v.calcification_rate = p.calcification_rate;

        v.relative_rupture_volume = p.relative_rupture_volume;
        v.rupture_volume = v.total * v.relative_rupture_volume;
    }
}