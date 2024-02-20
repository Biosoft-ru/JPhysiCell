package ru.biosoft.physicell.sample_projects.pred_prey_farmer;

import java.util.List;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

/* predator functions */
public class PredatorPhenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        CellDefinition pFarmerDef = CellDefinition.getCellDefinition( "farmer" );
        CellDefinition pPreyDef = CellDefinition.getCellDefinition( "prey" );
        CellDefinition pPredDef = CellDefinition.getCellDefinition( "predator" );

        // hunting 

        double max_detection_distance = 2;

        // see who is nearby 

        List<Cell> nearby = PredPreyFarmer.get_possible_neighbors( pCell );

        for( int i = 0; i < nearby.size(); i++ )
        {
            Cell pC = nearby.get( i );
            // is it prey ? 

            if( pC.type == pPreyDef.type )
            {
                boolean eat_it = true;
                // in range? 
                double[] displacement = VectorUtil.newDiff( pC.position, pCell.position );
                //			displacement -= pCell.position; 
                double distance = VectorUtil.norm( displacement );
                if( distance > pCell.phenotype.geometry.radius + pC.phenotype.geometry.radius + max_detection_distance )
                {
                    eat_it = false;
                }

                // am I hungry? 

                if( eat_it == true )
                {
                    // eat it! 
                    pCell.ingestCell( pC );

                    // increase energy 
                    pCell.custom_data.set( "energy", pCell.custom_data.get( "energy" ) + 100 );
                    //                pCell.custom_data["energy"] += 100; 	
                }
            }
        }

        // update energy 

        double decay_rate = 0.00025;
        pCell.custom_data.set( "energy", pCell.custom_data.get( "energy" ) / ( 1.0 + dt * decay_rate ) );
        //    pCell.custom_data["energy"] /= ( 1.0 + dt * decay_rate );

        // low energy kills

        // death based on food
        int nNecrosis = phenotype.death.findDeathModelIndex( "necrosis" );

        if( pCell.custom_data.get( "energy" ) < 0.1 )
        {
            pCell.startDeath( nNecrosis );
            pCell.functions.updatePhenotype = null;
            return;
        }
    }
}