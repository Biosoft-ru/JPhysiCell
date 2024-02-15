package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.Volume;

public class standard_volume_update_function implements CellFunctions.volume_update_function
{
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        Volume v = phenotype.volume;
        v.fluid += dt * v.fluid_change_rate * ( v.target_fluid_fraction * v.total - v.fluid );
        if( v.fluid < 0.0 )
            v.fluid = 0.0;

        v.nuclear_fluid = ( v.nuclear / ( v.total + 1e-16 ) ) * ( v.fluid );
        v.cytoplasmic_fluid = v.fluid - v.nuclear_fluid;

        v.nuclear_solid += dt * v.nuclear_biomass_change_rate * ( v.target_solid_nuclear - v.nuclear_solid );
        if( v.nuclear_solid < 0.0 )
            v.nuclear_solid = 0.0;

        v.target_solid_cytoplasmic = v.target_cytoplasmic_to_nuclear_ratio * v.target_solid_nuclear;// phenotype.volume.cytoplasmic_to_nuclear_fraction * 

        v.cytoplasmic_solid += dt * v.cytoplasmic_biomass_change_rate * ( v.target_solid_cytoplasmic - v.cytoplasmic_solid );
        if( v.cytoplasmic_solid < 0.0 )
            v.cytoplasmic_solid = 0.0;

        v.solid = v.nuclear_solid + v.cytoplasmic_solid;

        v.nuclear = v.nuclear_solid + v.nuclear_fluid;
        v.cytoplasmic = v.cytoplasmic_solid + v.cytoplasmic_fluid;

        v.calcified_fraction += dt * v.calcification_rate * ( 1 - v.calcified_fraction );

        v.total = v.cytoplasmic + v.nuclear;

        v.fluid_fraction = v.fluid / ( 1e-16 + v.total );

        phenotype.geometry.update( pCell, phenotype, dt );
    }
}