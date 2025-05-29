package ru.biosoft.physicell.covid;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.ui.AgentColorer;
import ru.biosoft.physicell.xml.ModelReader;

public class CovidColorer implements AgentColorer
{

    @Override
    public Color[] findColors(Cell pCell)
    {
        Model model = pCell.getModel();
        int lung_epithelial_type = model.getCellDefinition( "lung epithelium" ).type;

        int CD8_Tcell_type = model.getCellDefinition( "CD8 Tcell" ).type;
        int Macrophage_type = model.getCellDefinition( "macrophage" ).type;
        int Neutrophil_type = model.getCellDefinition( "neutrophil" ).type;
        int DC_type = model.getCellDefinition( "DC" ).type;
        int CD4_Tcell_type = model.getCellDefinition( "CD4 Tcell" ).type;
        int fibroblast_type = model.getCellDefinition( "fibroblast" ).type;
        int res_type = model.getCellDefinition( "residual" ).type;

        // start with white 
        Color[] output = new Color[] {Color.white, Color.black, Color.white, Color.white};
        //            std::vector<std::string> output = {"white", "black", "white" , "white" };   
        // false_cell_coloring_cytometry(pCell); 

        if( pCell.phenotype.death.dead == true )
        {
            if( pCell.type != lung_epithelial_type )
            {
                output[0] = ModelReader.readColor( model.getParameterString( "apoptotic_immune_color" ) );
                output[2] = output[0];
                output[3] = output[0];
                return output;
            }

            output[0] = ModelReader.readColor( model.getParameterString( "apoptotic_epithelium_color" ) );
            output[2] = output[0];
            output[3] = output[0];
            return output;
        }

        if( pCell.phenotype.death.dead == false && pCell.type == lung_epithelial_type )
        {
            // color by virion 
            output = epithelium_coloring_function( pCell );
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

            // (Adrianne) added colours to show when macrophages are exhausted and when they are hyperactivated
            if( pCell.phenotype.volume.total > pCell.customData.get( "threshold_macrophage_volume" ) )// macrophage exhausted
            {
                color = ModelReader.readColor( model.getParameterString( "exhausted_macrophage_color" ) );
            }
            else if( pCell.customData.get( "ability_to_phagocytose_infected_cell" ) == 1 )// macrophage has been activated to kill infected cells by T cell
            {
                color = ModelReader.readColor( model.getParameterString( "hyperactivated_macrophage_color" ) );
            }        

            output[0] = color;
            output[2] = output[0];
            output[3] = output[0];
            return output;
        }

        if( pCell.phenotype.death.dead == false && pCell.type == Neutrophil_type )
        {
            output[0] = ModelReader.readColor( model.getParameterString( "Neutrophil_color" ) );
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

        if( pCell.phenotype.death.dead == false && pCell.type == fibroblast_type )
        {
            output[0] = ModelReader.readColor( model.getParameterString( "fibroblast_color" ) );
            output[2] = output[0];
            output[3] = output[0];
            return output;
        }

        if( pCell.phenotype.death.dead == false && pCell.type == res_type )
        {
            output[0] = Color.black;
            output[2] = output[0];
            output[3] = output[0];
            return output;
        }

        return output;
    }

    Color[] epithelium_coloring_function(Cell pCell)
    {
        Color[] output = new Color[] {Color.black, Color.black, Color.black, Color.black};
        //            std::vector<std::string> output( 4, "black" ); /

        // static int color_index = cell_defaults.custom_data.find_variable_index( "assembled virion" ); 
        int color_index = pCell.customData.findVariableIndex( pCell.getModel().getParameterString( "color_variable" ) );
        int nV = pCell.customData.findVariableIndex( "virion" );

        // color by assembled virion 

        if( pCell.phenotype.death.dead == false )
        {
            // find fraction of max viral load 
            double v = pCell.customData.get( color_index );

            double interpolation = 0;
            if( v < 1 )
            {
                interpolation = 0;
            }
            if( v >= 1.0 && v < 10 )
            {
                interpolation = 0.25;
            }
            if( v >= 10.0 && v < 100 )
            {
                interpolation = 0.5;
            }
            if( v >= 100.0 && v < 1000 )
            {
                interpolation = 0.75;
            }
            if( v >= 1000.0 )
            {
                interpolation = 1.0;
            }

            int red = (int)Math.floor( 255.0 * interpolation );
            int green = red;
            int blue = 255 - red;

            Color c = new Color( red, green, blue );
            //                char color [1024]; 
            //                sprintf( color, "rgb(%u,%u,%u)" , red,green,blue ); 

            output[0] = c;
            output[2] = c;
            output[3] = c;
        }

        return output;
    }


}