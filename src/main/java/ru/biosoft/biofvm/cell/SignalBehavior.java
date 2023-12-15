package ru.biosoft.biofvm.cell;

import java.util.HashMap;
import java.util.Map;

import ru.biosoft.biofvm.Microenvironment;
import ru.biosoft.biofvm.VectorUtil;

public class SignalBehavior
{

    public static Map<String, Integer> signal_to_int = new HashMap<>();
    public static double[] signal_scales = new double[0];
    public static double get_single_signal( Cell pCell, int index )
    {
        Microenvironment microenvironment = Microenvironment.get_default_microenvironment();
        int m = microenvironment.number_of_densities(); 
        int n = Cell.cell_definition_indices_by_name.size(); 

        double out = 0.0; 
        if( index < 0 )
        { 
//            std::cout<< "Why would you ask for array[-1]? Why? WHY???? That's it, I quit." << std::endl; 
            return -9e9; 
        }

        // first m entries: extracellular concentration 
        int start_substrate_ind = find_signal_index( microenvironment.density_names[0] ); 
        if( start_substrate_ind <= index && index < start_substrate_ind + m )
        {
            out = pCell.nearest_density_vector()[index-start_substrate_ind];
            out /= signal_scales[index]; 
            return out; 
        }

        // second m entries: intracellular concentration 
        int start_int_substrate_ind = find_signal_index( "intracellular " + microenvironment.density_names[0] ); 
        if( start_int_substrate_ind <= index && index < start_int_substrate_ind + m )
        {
            out = pCell.phenotype.molecular.internalized_total_substrates[index-start_int_substrate_ind]; 
            out /= pCell.phenotype.volume.total;
            out /= signal_scales[index]; 
            return out; 
        }

        // next m entries: gradients 
        int start_substrate_grad_ind = find_signal_index( microenvironment.density_names[0] + " gradient"); 
        if( start_substrate_grad_ind <= index && index < start_substrate_grad_ind + m )
        {
            out = VectorUtil.norm( pCell.nearest_gradient( index - start_substrate_grad_ind ) );
            out /= signal_scales[index]; 
            return out; 
        }

        // mechanical pressure 
        int pressure_ind = find_signal_index( "pressure" ); 
        if( index == pressure_ind )
        {
            out = pCell.state.simple_pressure;
            out /= signal_scales[index]; 
            return out; 
        }

        // cell volume  
        int volume_ind = find_signal_index( "volume"); 
        if( index == volume_ind )
        {
            out = pCell.phenotype.volume.total; 
            out /= signal_scales[index]; 
            return out; 
        }

        // physical contact with cells (of each type) 
        // individual contact signals are a bit costly 
        int contact_ind = find_signal_index( "contact with " + Cell.cell_definitions_by_type.get(0).name ); 
        if( contact_ind <= index && index < contact_ind + n+2 )
        {
            //            std::vector<int> counts( n , 0 );
            int[] counts = new int[n];
            // process all neighbors 
            int dead_cells = 0; 
            int live_cells = 0; 
            for( Cell pC : pCell.state.neighbors )//int i=0; i < pCell.state.neighbors.size(); i++ )
            {
                //                Cell pC = pCell.state.neighbors[i];
                if( pC.phenotype.death.dead == true )
                { dead_cells++; } 
                else
                { live_cells++; } 
                int nCT = Cell.cell_definition_indices_by_type.get( pC.type );
                counts[nCT] += 1; 
            }

            if( index < contact_ind + n )
            {
                out = counts[index-contact_ind]; 
                out /= signal_scales[index]; 
                return out; 
            }

            int live_contact_ind = find_signal_index( "contact with live cell"); 
            if( index == live_contact_ind )
            {
                out = live_cells; 
                out /= signal_scales[index]; 
                return out; 
            }

            int dead_contact_ind = find_signal_index( "contact with dead cell"); 
            // index == dead_contact_ind
            out = dead_cells; 
            out /= signal_scales[index]; 
            return out; 
        }

        // contact with BM 
        int BM_contact_ind = find_signal_index( "contact with basement membrane"); 
        if( index == BM_contact_ind )
        { 
            //            out = (double) pCell.state.contact_with_basement_membrane; 
            out = pCell.state.contact_with_basement_membrane ? 1 : 0;
            out /= signal_scales[index]; 
            return out; 
        } 

        // damage
        int damage_ind = find_signal_index( "damage"); 
        if( index == damage_ind )
        {
            out = pCell.state.damage; 
            out /= signal_scales[index]; 
            return out; 
        } 

        // live / dead state 
        int dead_ind = find_signal_index( "dead" ); 
        if( index == dead_ind )
        {
            //            out = (double) pCell.phenotype.death.dead;  
            out = pCell.phenotype.death.dead ? 1 : 0;
            out /= signal_scales[index]; 
            return out; 
        } 

        // integrated total attack time 
        int tot_attack_ind = find_signal_index( "total attack time"); 
        if( index == tot_attack_ind )
        {
            out = pCell.state.total_attack_time;     
            out /= signal_scales[index]; 
            return out; 
        } 

        // time 
        int time_ind = find_signal_index( "time" ); 
        if( index == time_ind )
        {
            //            out = PhysiCell_globals.current_time;      
            //            out /= signal_scales[index]; //TODO: later
            //            return out; 
        } 

        // custom signals 
        int first_custom_ind = find_signal_index( "custom 0" ); 
        int max_custom_ind = first_custom_ind + pCell.custom_data.variables.size();
        if( first_custom_ind > -1 && index >= first_custom_ind && index < max_custom_ind )
        {
            out = pCell.custom_data.variables.get( index - first_custom_ind ).value;
            out /= signal_scales[index];
            return out; 
        }

        int apoptotic_ind = find_signal_index( "apoptotic" ); 
        if( index == apoptotic_ind )
        {
            if( pCell.phenotype.cycle.current_phase().code == PhysiCellConstants.apoptotic )
            { return 1; }
            else
            { return 0; }
        }

        int necrotic_ind = find_signal_index( "necrotic" ); 
        if( index == necrotic_ind )
        {
            if( pCell.phenotype.cycle.current_phase().code == PhysiCellConstants.necrotic_swelling
                    || pCell.phenotype.cycle.current_phase().code == PhysiCellConstants.necrotic_lysed
                    || pCell.phenotype.cycle.current_phase().code == PhysiCellConstants.necrotic )
            { return 1; }
            else
            { return 0; }
        }

    /*
        static int start_immunogenicity_ind = find_signal_index( "immunogenicity to " + cell_definitions_by_type[0].name ); 
        static int max_immunogenicity_ind = start_immunogenicity_ind+n; 
        if( start_immunogenicity_ind > -1 && index >= start_immunogenicity_ind && index < max_immunogenicity_ind )
        {
            int j = index - start_immunogenicity_ind; 
            out = pCell.phenotype.cell_interactions.immunogenicities[j]; 
            out /= signal_scales[index];
            return out; 
        }
    */


        // unknown after here !

        System.out.println( "Warning: Requested unknown signal number " + index + "!" );
        System.out.println("         Returning 0.0, but you should fix this!");
//        std::cout << "Warning: Requested unknown signal number " << index << "!" << std::endl
//                  << "         Returning 0.0, but you should fix this!" << std::endl << std::endl; 
        return 0.0; 
    }

    static int find_signal_index(String signal_name)
    {
        //        auto search = signal_to_int.find( signal_name );
        // safety first! 
        //        if( search != signal_to_int.end() )
        //        {
        //            return search -> second;
        //        }
        if( signal_to_int.containsKey( signal_name ) )
            return signal_to_int.get( signal_name );

        System.out.println( "having trouble finding " + signal_name );
        //        std::cout << "having trouble finding " << signal_name << std::endl; 

        return -1;
    }

    public static double get_single_signal(Cell pCell, String name)
    {
        return get_single_signal( pCell, find_signal_index( name ) );
    }
}
