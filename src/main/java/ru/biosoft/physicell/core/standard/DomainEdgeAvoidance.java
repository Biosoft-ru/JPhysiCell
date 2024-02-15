package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.CellMembraneInteractions;

public class DomainEdgeAvoidance implements CellMembraneInteractions
{
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        if( pCell.functions.calculate_distance_to_membrane == null )
        {
            pCell.functions.calculate_distance_to_membrane = new DomainEdgeDistance();
        }
        phenotype.mechanics.cellBMRepulsionStrength = 100;

        double max_interactive_distance = phenotype.mechanics.relMaxAdhesionDistance * phenotype.geometry.radius;
        double distance = pCell.functions.calculate_distance_to_membrane.execute( pCell, phenotype, dt );
        //Note that the distance_to_membrane function must set displacement values (as a normal vector)

        // Repulsion from basement membrane
        double temp_r = 0;
        if( distance < phenotype.geometry.radius )
        {
            temp_r = ( 1 - distance / phenotype.geometry.radius );
            temp_r *= temp_r;
            temp_r *= phenotype.mechanics.cellBMRepulsionStrength;
        }
        if( Math.abs( temp_r ) < 1e-16 )
            return;

        VectorUtil.axpy( ( pCell.velocity ), temp_r, pCell.displacement );
    }
}