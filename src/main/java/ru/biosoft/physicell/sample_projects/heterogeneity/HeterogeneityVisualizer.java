package ru.biosoft.physicell.sample_projects.heterogeneity;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.ui.AgentVisualizer;

public class HeterogeneityVisualizer extends AgentVisualizer
{
    double p_min;
    double p_max;

    public HeterogeneityVisualizer(Model model)
    {
        p_min = model.getParameterDouble( "oncoprotein_min" );
        p_max = model.getParameterDouble( "oncoprotein_max" );
    }

    @Override
    public Color findColor(Cell pCell)
    {
        double p = SignalBehavior.getSingleSignal( pCell, "custom:oncoprotein" );

        // immune are black
        Color output = Color.black;

        if( pCell.type == 1 )
        {
            return output;
        }

        // live cells are green, but shaded by oncoprotein value 
        if( pCell.phenotype.death.dead == false )
        {
            int oncoprotein = (int)Math.round( ( 1.0 / ( p_max - p_min ) ) * ( p - p_min ) * 255.0 );
            output = new Color( (int) ( oncoprotein / 1 ), (int) ( oncoprotein / 1 ), (int) ( ( 255 - oncoprotein ) / 1 ) );
            //                output = new Color( (int) ( oncoprotein / p_max ), (int) ( oncoprotein / p_max ), (int) ( ( 255 - oncoprotein ) / p_max ) );
            //		char szTempString [128];
            //		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
            //		output[0].assign( szTempString );
            //		output[1].assign( szTempString );
            //
            //		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/p_max) , (int)round(output[0][1]/p_max) , (int)round(output[0][2]/p_max) );
            //		output[2].assign( szTempString );
            return output;
        }

        // if not, dead colors 
        if( SignalBehavior.getSingleSignal( pCell, "apoptotic" ) > 0.5 )
        {
            output = new Color( 125, 0, 0 );
            //		output[0] = "rgb(255,0,0)";
            //		output[2] = "rgb(125,0,0)";
        }

        // Necrotic - Brown
        if( SignalBehavior.getSingleSignal( pCell, "necrotic" ) > 0.5 )
        {
            output = new Color( 139, 69, 19 );
            //		output[0] = "rgb(250,138,38)";
            //		output[2] = "rgb(139,69,19)";
        }
        return output;
    }
}