package ru.biosoft.physicell.sample_projects.pred_prey_farmer;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.ui.AgentVisualizer2;

public class PPFVisualizer extends AgentVisualizer2
{
    @Override
    public Color[] findColors(Cell pCell)
    {
        CellDefinition pFarmerDef = CellDefinition.getCellDefinition( "farmer" );
        CellDefinition pPreyDef = CellDefinition.getCellDefinition( "prey" );
        CellDefinition pPredDef = CellDefinition.getCellDefinition( "predator" );

        if( pCell.type == pFarmerDef.type )
            return new Color[] {Color.gray};

        if( pCell.type == pPreyDef.type )
            return new Color[] {Color.blue};

        if( pCell.type == pPredDef.type )
            return new Color[] {Color.orange};
        return super.findColors( pCell );
    }
}