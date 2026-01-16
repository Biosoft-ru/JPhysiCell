package ru.biosoft.physicell.micromets;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.ui.AgentColorer;
import ru.biosoft.physicell.xml.ModelReader;

public class MicrometsColorer implements AgentColorer
{

    @Override
    public Color[] findColors(Cell pCell)
    {
        Model model = pCell.getModel();
        int lung_epithelial_type = model.getCellDefinition( "lung cell" ).type;
        int cancer_type = model.getCellDefinition( "cancer cell" ).type;

        int CD8_Tcell_type = model.getCellDefinition( "CD8 Tcell" ).type;
        int Macrophage_type = model.getCellDefinition( "macrophage" ).type;
        int DC_type = model.getCellDefinition( "DC" ).type;
        int CD4_Tcell_type = model.getCellDefinition( "CD4 Tcell" ).type;

        // start with white

        Color[] output = {Color.white, Color.black, Color.white, Color.white};
        // false_cell_coloring_cytometry(pCell);

        if( pCell.phenotype.death.dead == true )
        {
            if( pCell.type != lung_epithelial_type && pCell.type != cancer_type )
            {
                output[0] = ModelReader.readColor( model.getParameterString( "apoptotic_immune_color" ) );
                output[2] = output[0];
                output[3] = output[0];
                return output;
            }
            else
            {
                if( pCell.type == lung_epithelial_type )
                {
                    output[0] = ModelReader.readColor( model.getParameterString( "apoptotic_epithelium_color" ) );
                    output[2] = output[0];
                    output[3] = output[0];
                    return output;
                }
                else
                {
                    output[0] = ModelReader.readColor( model.getParameterString( "apoptotic_cancer_color" ) );
                    output[2] = output[0];
                    output[3] = output[0];
                }
            }
        }

        if( pCell.phenotype.death.dead == false && pCell.type == lung_epithelial_type )
        {
            output[0] = ModelReader.readColor( "rgb(0,0,255)" );
            output[2] = ModelReader.readColor( "rgb(0,0,255)" );
            output[3] = ModelReader.readColor( "rgb(0,0,255)" );
            // output = epithelium_coloring_function(pCell);
            return output;
        }

        if( pCell.phenotype.death.dead == false && pCell.type == cancer_type )
        {
            output[0] = ModelReader.readColor( model.getParameterString( "cancer_color" ) );
            output[2] = output[0];
            output[3] = output[0];
            return output;
        }

        if( pCell.phenotype.death.dead == false && pCell.type == CD8_Tcell_type )
        {
            output[0] = ModelReader.readColor( model.getParameterString( "CD8_Tcell_color" ) );
            output[2] = output[0];
            output[3] = output[0];
            return output;
        }

        // (Adrianne) adding CD4 T cell colouring
        if( pCell.phenotype.death.dead == false && pCell.type == CD4_Tcell_type )
        {
            output[0] = ModelReader.readColor( model.getParameterString( "CD4_Tcell_color" ) );
            output[2] = output[0];
            output[3] = output[0];
            return output;
        }

        if( pCell.phenotype.death.dead == false && pCell.type == Macrophage_type )
        {
            Color color = ModelReader.readColor( model.getParameterString( "Macrophage_color" ) );
            if( pCell.customData.get( "activated_immune_cell" ) > 0.5 )
            {
                color = ModelReader.readColor( model.getParameterString( "activated_macrophage_color" ) );
            }

            if( pCell.phenotype.volume.total > pCell.customData.get( "threshold_macrophage_volume" ) )// macrophage exhausted
            {
                color = ModelReader.readColor( model.getParameterString( "exhausted_macrophage_color" ) );
            }
            else if( pCell.customData.get( "ability_to_phagocytose_cancer_cell" ) == 1 )// macrophage has been activated to kill cancer cells by T cell
            {
                color = ModelReader.readColor( model.getParameterString( "hyperactivated_macrophage_color" ) );
            }

            output[0] = color;
            output[2] = output[0];
            output[3] = output[0];
            return output;
        }

        //(Adrianne) adding colour for DCs
        if( pCell.phenotype.death.dead == false && pCell.type == DC_type )
        {
            Color color = ModelReader.readColor( model.getParameterString( "DC_color" ) );
            if( pCell.customData.get( "activated_immune_cell" ) > 0.5 )
            {
                color = ModelReader.readColor( model.getParameterString( "activated_DC_color" ) );
            }

            output[0] = color;
            output[2] = output[0];
            output[3] = output[0];
            return output;
        }

        return output;
    }
}