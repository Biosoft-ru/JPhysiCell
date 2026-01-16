package ru.biosoft.physicell.micromets;

import java.util.Set;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.RandomGenerator;

public class Util
{

    public static boolean check_for_out_of_bounds(Cell pC, double tolerance) // return true if out of bounds, within a tolerance
    {
        Microenvironment microenvironment = pC.getMicroenvironment();
        double Xmin = microenvironment.mesh.boundingBox[0];
        double Ymin = microenvironment.mesh.boundingBox[1];
        double Zmin = microenvironment.mesh.boundingBox[2];

        double Xmax = microenvironment.mesh.boundingBox[3];
        double Ymax = microenvironment.mesh.boundingBox[4];
        double Zmax = microenvironment.mesh.boundingBox[5];

        boolean two_dimensions = pC.getMicroenvironment().getOptions().simulate2D;

        boolean setup_done = false;
        if( pC.getMicroenvironment().getOptions().simulate2D && setup_done == false )
        {
            Zmin = 0.0;
            Zmax = 0.0;
            setup_done = true;
        }

        if( pC.position[0] < Xmin + tolerance )
        {
            return true;
        }
        if( pC.position[0] > Xmax - tolerance )
        {
            return true;
        }

        if( pC.position[1] < Ymin + tolerance )
        {
            return true;
        }
        if( pC.position[1] > Ymax - tolerance )
        {
            return true;
        }

        if( two_dimensions )
        {
            return false;
        }

        if( pC.position[2] < Zmin + tolerance )
        {
            return true;
        }
        if( pC.position[2] > Zmax - tolerance )
        {
            return true;
        }

        return false;
    }
    
    public static Cell immune_cell_check_neighbors_for_attachment(Cell pAttacker, double dt)
    {
        int CD8_Tcell_type = pAttacker.getModel().getCellDefinition( "CD8 Tcell" ).type;
        Set<Cell> nearby = pAttacker.cells_in_my_container();
        //        int i = 0;
        //        while( i < nearby.size() )
        //        {
        for( Cell nearbyCell : nearby )
        {
            //             don't try to kill yourself
            if( nearbyCell != pAttacker )
            {
                // CD8 T cell attaches only to one cell
                if( attempt_immune_cell_attachment( pAttacker, nearbyCell, dt ) && pAttacker.type == CD8_Tcell_type )
                {
                    return nearbyCell;
                }
            }
            //            i++;
        }

        // Cell has more than one attached cell
        if( pAttacker.state.numberAttachedCells() > 0 )
            return pAttacker.state.attachedCells.iterator().next();

        return null;
    }

    public static boolean attempt_immune_cell_attachment(Cell pAttacker, Cell pTarget, double dt)
    {
        Model model = pAttacker.getModel();
        RandomGenerator rng = model.getRNG();
        int CD8_Tcell_type = model.getCellDefinition( "CD8 Tcell" ).type;
        int DC_type = model.getCellDefinition( "DC" ).type;
        int cancer_type = model.getCellDefinition( "cancer cell" ).type;
        int lung_type = model.getCellDefinition( "lung cell" ).type;

        // if the target is not cancer cell, give up for CD8 attack
        if( pTarget.type != cancer_type && pAttacker.type == CD8_Tcell_type )
        {
            return false;
        }

        // If the target is not cancer cell, give up for DC attack (Just cancer cells attach DCs)
        //if ( pTarget.type != cancer_type  && pAttacker.type == DC_type) //

        // If the target is not cancer or lung cell, give up for DC attack (Just cancer or lung cells attach DCs)
        if( pTarget.type != cancer_type && pTarget.type != lung_type && pAttacker.type == DC_type )
        {
            return false;
        }

        // if the target cell is dead, give up for CD8 attack
        if( pTarget.phenotype.death.dead == true && pAttacker.type == CD8_Tcell_type )
        {
            return false;
        }

        // if the target cell is not dead, give up for DC attach (cancer dead cell presenting antigen)
        if( pTarget.phenotype.death.dead != true && pAttacker.type == DC_type )
        {
            return false;
        }

        // if the target cell is too far away, give up

        double[] displacement = VectorUtil.newDiff( pTarget.position, pAttacker.position );
        double distance_scale = VectorUtil.norm( displacement );

        // better: use mechanics constants
        if( distance_scale > pAttacker.customData.get( "max_attachment_distance" ) )
        {
            return false;
        }

        // now, get the attachment probability
        double attachment_probability = pAttacker.customData.get( "cell_attachment_rate" ) * dt;

        // don't need to cap it at 1.00: if prob > 100%,
        // then this statement always evaluates as true,
        // just the same as capping probability at 100%
        if( rng.UniformRandom() <= attachment_probability )
        {
            Cell.attachcCells( pAttacker, pTarget );
            return true;
        }

        return false;
    }

}
