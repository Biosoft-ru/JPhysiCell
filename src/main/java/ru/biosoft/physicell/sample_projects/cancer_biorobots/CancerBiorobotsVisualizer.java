package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.ui.AgentVisualizer;

public class CancerBiorobotsVisualizer extends AgentVisualizer
{
    @Override
    public Color findBorderColor(Cell cell)
    {
        return Color.black;
    }

    @Override
    public Color[] findColors(Cell cell)
    {
        double damage = SignalBehavior.getSingleSignal( cell, "damage" );
        double maxDamage = 1.0 * SignalBehavior.getSingleSignal( cell, "custom:damage_rate" )
                / ( 1e-16 + SignalBehavior.getSingleSignal( cell, "custom:repair_rate" ) );

        CellDefinition pCD_cargo = cell.getModel().getCellDefinition( "cargo cell" );
        CellDefinition pCD_worker = cell.getModel().getCellDefinition( "worker cell" );

        if( cell.type == pCD_cargo.type )
            return new Color[] {Color.blue, Color.blue};

        if( cell.type == pCD_worker.type )
            return new Color[] {Color.red, Color.red};

        Color[] output = new Color[] {Color.black, Color.black, Color.black, Color.black};

        if( SignalBehavior.getSingleSignal( cell, "apoptotic" ) > 0.5 ) // Apoptotic - cyan
        {
            output[0] = Color.cyan;
            output[2] = Color.cyan.darker();
            return output;
        }

        if( SignalBehavior.getSingleSignal( cell, "necrotic" ) > 0.5 )
        {
            output[0] = new Color( 250, 138, 38 );
            output[2] = new Color( 139, 69, 19 );
            return output;
        }
        // live tumor -- shade by level of damage  if live: color by damage 
        if( SignalBehavior.getSingleSignal( cell, "dead" ) < 0.5 )
        {
            int damageInt = (int)Math.round( damage * 255.0 / maxDamage );
            Color c = new Color( damageInt, 255 - damageInt, damageInt );
            Color c2 = new Color( damageInt / 4, ( 255 - damageInt ) / 4, damageInt / 4 );
            return new Color[] {c, c, c2, c2};
        }
        return output;
    }
}