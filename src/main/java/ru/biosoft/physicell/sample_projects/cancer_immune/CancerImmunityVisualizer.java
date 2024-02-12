package ru.biosoft.physicell.sample_projects.cancer_immune;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.ui.AgentVisualizer;

public class CancerImmunityVisualizer extends AgentVisualizer
{
    @Override
    public Color findBorderColor(Cell cell)
    {
        return Color.black;
    }

    @Override
    public Color findColor(Cell cell)
    {
        int oncoprotein_i = cell.custom_data.find_variable_index( "oncoprotein" );
        Color result = Color.black; // immune are black

        if( cell.type == 1 )
        {
            return new Color( 50, 205, 50 );
        }

        // if I'm under attack, color me 
        if( cell.state.attachedCells.size() > 0 )
        {
            return new Color( 128, 0, 128 );
        }
        // live cells are green, but shaded by oncoprotein value 

        if( cell.phenotype.death.dead == false )
        {
            int oncoprotein = (int)Math.round( 0.5 * cell.custom_data.get( oncoprotein_i ) * 255.0 );
            return new Color( oncoprotein, oncoprotein, 255 - oncoprotein );
            //                        char szTempString [128];
            //                        sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
            //                        output[0].assign( szTempString );
            //                        output[1].assign( szTempString );

            //                        sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
            //                        output[2].assign( szTempString );

            //                        return output; 
        }

        // if not, dead colors 

        if( cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.apoptotic ) // Apoptotic - Red
        {
            return new Color( 255, 0, 0 );
            //            output[0] = "rgb(255,0,0)";
            //            output[2] = "rgb(125,0,0)";
        }

        //        // Necrotic - Brown
        if( cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                || cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                || cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
        {
            return new Color( 250, 138, 138 );
            //                        output[0] = "rgb(250,138,38)";
            //                        output[2] = "rgb(139,69,19)";
        }
        //        
        return result;
    }
}