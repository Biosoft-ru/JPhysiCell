package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.MembraneInteractions;
import ru.biosoft.physicell.core.Phenotype;

public class DomainEdgeAvoidance extends MembraneInteractions
{
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        if( pCell.functions.membraneDistanceCalculator == null )
        {
            pCell.functions.membraneDistanceCalculator = new DomainEdgeDistance();
        }
        phenotype.mechanics.cellBMRepulsionStrength = 100;

        double distance = pCell.functions.membraneDistanceCalculator.execute( pCell, phenotype, dt );
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


    @Override
    public String getName()
    {
        return "Avoid domain edge";
    }
}