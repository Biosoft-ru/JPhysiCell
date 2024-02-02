package ru.biosoft.physicell.core;

import java.util.HashMap;
import java.util.Map;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.CustomCellData.Variable;

public class SignalBehavior
{
    public static Map<String, Integer> behavior_to_int = new HashMap<>();
    public static Map<Integer, String> int_to_behavior = new HashMap<>();
    static Map<Integer, String> int_to_signal = new HashMap<>();
    public static Map<String, Integer> signal_to_int = new HashMap<>();
    public static double[] signal_scales = new double[0];

    private static boolean setupDone = false;
    public static double get_single_signal(Cell pCell, int index)
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();//Microenvironment.get_default_microenvironment();
        int m = microenvironment.number_of_densities();
        int n = CellDefinition.getDefinitionsCount();

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
            out = pCell.nearest_density_vector()[index - start_substrate_ind];
            out /= signal_scales[index];
            return out;
        }

        // second m entries: intracellular concentration 
        int start_int_substrate_ind = find_signal_index( "intracellular " + microenvironment.density_names[0] );
        if( start_int_substrate_ind <= index && index < start_int_substrate_ind + m )
        {
            out = pCell.phenotype.molecular.internalized_total_substrates[index - start_int_substrate_ind];
            out /= pCell.phenotype.volume.total;
            out /= signal_scales[index];
            return out;
        }

        // next m entries: gradients 
        int start_substrate_grad_ind = find_signal_index( microenvironment.density_names[0] + " gradient" );
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
            out = pCell.state.simplePressure;
            out /= signal_scales[index];
            return out;
        }

        // cell volume  
        int volume_ind = find_signal_index( "volume" );
        if( index == volume_ind )
        {
            out = pCell.phenotype.volume.total;
            out /= signal_scales[index];
            return out;
        }

        // physical contact with cells (of each type) 
        // individual contact signals are a bit costly 
        int contact_ind = find_signal_index( "contact with " + CellDefinition.getCellDefinition( 0 ).name );
        if( contact_ind <= index && index < contact_ind + n + 2 )
        {
            int[] counts = new int[n];
            // process all neighbors 
            int dead_cells = 0;
            int live_cells = 0;
            for( Cell pC : pCell.state.neighbors )
            {
                //                Cell pC = pCell.state.neighbors[i];
                if( pC.phenotype.death.dead == true )
                {
                    dead_cells++;
                }
                else
                {
                    live_cells++;
                }
                //                int nCT = CellDefinition.getCellDefinition( pC.type );// Cell.cell_definition_indices_by_type.get( pC.type );
                counts[pC.type] += 1;
            }

            if( index < contact_ind + n )
            {
                out = counts[index - contact_ind];
                out /= signal_scales[index];
                return out;
            }

            int live_contact_ind = find_signal_index( "contact with live cell" );
            if( index == live_contact_ind )
            {
                out = live_cells;
                out /= signal_scales[index];
                return out;
            }

            int dead_contact_ind = find_signal_index( "contact with dead cell" );
            // index == dead_contact_ind TODO: check
            out = dead_cells;
            out /= signal_scales[index];
            return out;
        }

        // contact with BM 
        int BM_contact_ind = find_signal_index( "contact with basement membrane" );
        if( index == BM_contact_ind )
        {
            //            out = (double) pCell.state.contact_with_basement_membrane; 
            out = pCell.state.contactWithBasementMembrane ? 1 : 0;
            out /= signal_scales[index];
            return out;
        }

        // damage
        int damage_ind = find_signal_index( "damage" );
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
            out = pCell.phenotype.death.dead ? 1 : 0;
            out /= signal_scales[index];
            return out;
        }

        // integrated total attack time 
        int tot_attack_ind = find_signal_index( "total attack time" );
        if( index == tot_attack_ind )
        {
            out = pCell.state.totalAttackTime;
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
            if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.apoptotic )
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }

        int necrotic_ind = find_signal_index( "necrotic" );
        if( index == necrotic_ind )
        {
            if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                    || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                    || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
            {
                return 1;
            }
            else
            {
                return 0;
            }
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
        System.out.println( "         Returning 0.0, but you should fix this!" );
        return 0.0;
    }

    static int find_signal_index(String signal_name)
    {
        if( signal_to_int.containsKey( signal_name ) )
            return signal_to_int.get( signal_name );
        System.out.println( "having trouble finding " + signal_name );
        return -1;
    }

    public static double get_single_signal(Cell pCell, String name)
    {
        return get_single_signal( pCell, find_signal_index( name ) );
    }

    private static void register(String name, int index)
    {
        signal_to_int.put( name, index );// ] = m+i;
        int_to_signal.put( index, name );//] = name; 
    }

    private static void registerBehavior(String name, int index)
    {
        behavior_to_int.put( name, index );// ] = m+i;
        int_to_behavior.put( index, name );//] = name; 
    }

    public static void setup_signal_behavior_dictionaries(Microenvironment microenvironment)
    {
        // set key parameters on number of signals, etc. 
        // make registry of signals 
        // make registry of responses 
        if( setupDone == true )
            return;
        setupDone = true;

        int m = microenvironment.number_of_densities();

        signal_to_int.clear();
        int_to_signal.clear();

        // construct signals 

        // substrate densities 
        for( int i = 0; i < m; i++ )
        {
            register( microenvironment.density_names[i], i );
        }

        // internalized substrates 
        int map_index = m;
        for( int i = 0; i < m; i++ )
        {
            register( "intracellular " + microenvironment.density_names[i], m + i );
            signal_to_int.put( "internalized " + microenvironment.density_names[i], m + i );//[ name ] = m+i; 
        }

        // substrate gradients 
        map_index = 2 * m;
        for( int i = 0; i < m; i++ )
        {
            register( microenvironment.density_names[i] + " gradient", map_index );
            signal_to_int.put( "grad(" + microenvironment.density_names[i] + ")", map_index );//[ name ] = map_index;
            signal_to_int.put( "gradient of " + microenvironment.density_names[i], map_index );// ] = map_index;
            map_index++;
        }

        // mechanical pressure 
        register( "pressure", map_index );

        // total volume 
        map_index++;
        register( "volume", map_index );

        // contact with each cell type 
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            map_index++;
            signal_to_int.put( "contact with " + cd.name, map_index );
            register( "contact with cell type " + cd.type, map_index );
        }

        // contact with (any) live cell 
        map_index++;
        register( "contact with live cell", map_index );
        signal_to_int.put( "contact with live cells", map_index );

        // contact with dead cell 
        map_index++;
        register( "contact with dead cell", map_index );
        signal_to_int.put( "contact with dead cells", map_index );

        // contact with basement membrane 
        map_index++;
        register( "contact with basement membrane", map_index );
        signal_to_int.put( "contact with BM", map_index );

        // damage state 
        map_index++;
        register( "damage", map_index );

        // live / dead state 
        map_index++;
        register( "dead", map_index );
        signal_to_int.put( "is dead", map_index );

        // total attack time 
        map_index++;
        register( "total attack time", map_index );

        // current time
        map_index++;
        register( "time", map_index );
        signal_to_int.put( "current time", map_index );
        signal_to_int.put( "global time", map_index );

        int nc = 0;
        //        Set<String> customVars = new HashSet<>();
        //        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        //        {
        //            for( Variable var : cd.custom_data.variables )
        //                customVars.add( var.name );
        //        }
        CellDefinition def = CellDefinition.getCellDefinition( 0 );
        // custom signals
        for( Variable var : def.custom_data.variables )
        {
            String varName = var.name;
            map_index++;
            register( "custom:" + varName, map_index );
            signal_to_int.put( "custom: " + varName, map_index );
            signal_to_int.put( "custom " + nc, map_index );
            nc++;
        }

        map_index++;
        register( "apoptotic", map_index );
        signal_to_int.put( "is_apoptotic", map_index );

        map_index++;
        register( "necrotic", map_index );
        signal_to_int.put( "is_necrotic", map_index );

        /*
        // immunogenicity to each cell type 
        for( int i=0; i < n ; i++ )
        {
            map_index++; 
            Cell_Definition* pCD = cell_definitions_by_type[i]; 
            std::string temp =  "immunogenicity to " + pCD->name; 
            signal_to_int[temp] = map_index; 
            int_to_signal[map_index] = temp;        
                    // synonyms 
            std::string temp1 = "immunogenicity to cell type " + std::to_string( pCD->type ); 
            signal_to_int[temp1] = map_index; 
        }
        */

        /* add new signals above this line */
        behavior_to_int.clear();
        int_to_behavior.clear();

        // construct behaviors 
        String name;
        for( int i = 0; i < m; i++ )
        {
            map_index = i;
            name = microenvironment.density_names[i];

            // secretion rate 
            registerBehavior( name + " " + "secretion", map_index );

            // secretion target 
            map_index = m + i;
            registerBehavior( name + " " + "secretion target", map_index );

            // synonym 
            behavior_to_int.put( name + " " + "secretion saturation density", map_index );

            // uptake rate 
            map_index = 2 * m + i;
            registerBehavior( name + " " + "uptake", map_index );

            // net export rate 
            map_index = 3 * m + i;
            registerBehavior( name + " " + "export", map_index );
        }

        map_index = 4 * m;
        registerBehavior( "cycle entry", map_index );

        // synonym 
        behavior_to_int.put( "exit from cycle phase 0", map_index );

        // other cyle phases 
        for( int i = 1; i < 6; i++ )
        {
            map_index++;
            registerBehavior( "exit from cycle phase " + i, map_index );
        }

        map_index++;
        registerBehavior( "apoptosis", map_index );

        map_index++;
        registerBehavior( "necrosis", map_index );

        map_index++;
        registerBehavior( "migration speed", map_index );

        map_index++;
        registerBehavior( "migration bias", map_index );

        map_index++;
        registerBehavior( "migration persistence time", map_index );

        // chemotactic sensitivities 
        for( int i = 0; i < m; i++ )
        {
            map_index++;
            registerBehavior( "chemotactic response to " + microenvironment.density_names[i], map_index );
            behavior_to_int.put( "chemotactic sensitivity to " + microenvironment.density_names[i], map_index );
        }

        // cell-cell adhesion 
        map_index++;
        registerBehavior( "cell-cell adhesion", map_index );

        map_index++;
        registerBehavior( "cell-cell adhesion elastic constant", map_index );

        // cell adhesion affinities 
        // cell-type specific adhesion 
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            map_index++;
            registerBehavior( "adhesive affinity to " + cd.name, map_index );
            behavior_to_int.put( "adhesive affinity to cell type " + cd.type, map_index );
        }

        // max adhesion distance 
        map_index++;
        registerBehavior( "relative maximum adhesion distance", map_index );

        // cell-cell repulsion 
        map_index++;
        registerBehavior( "cell-cell repulsion", map_index );

        // cell-BM adhesion 
        map_index++;
        registerBehavior( "cell-BM adhesion", map_index );
        behavior_to_int.put( "cell-membrane adhesion", map_index );

        // cell-BM repulsion 
        map_index++;
        registerBehavior( "cell-BM repulsion", map_index );
        behavior_to_int.put( "cell-membrane repulsion", map_index );

        // phagocytosis of dead cell
        map_index++;
        registerBehavior( "phagocytose dead cell", map_index );
        behavior_to_int.put( "phagocytosis of dead cell", map_index );
        behavior_to_int.put( "phagocytosis of dead cells", map_index );

        // phagocytosis of each live cell type 
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            map_index++;
            registerBehavior( "phagocytose " + cd.name, map_index );
            behavior_to_int.put( "phagocytose cell type " + cd.type, map_index );
            behavior_to_int.put( "phagocytosis of " + cd.type, map_index );
        }

        // attack of each live cell type 
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            map_index++;
            registerBehavior( "attack " + cd.name, map_index );
            behavior_to_int.put( "attack cell type " + cd.type, map_index );
        }

        // fusion 
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            map_index++;
            registerBehavior( "fuse to " + cd.name, map_index );
            behavior_to_int.put( "fuse to cell type " + cd.type, map_index );
        }

        // transformation 
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            map_index++;
            registerBehavior( "transform to " + cd.name, map_index );
            behavior_to_int.put( "transform to cell type " + cd.type, map_index );
        }

        // custom behaviors
        nc = 0;
        for( Variable var : def.custom_data.variables )
        {
            String varName = var.name;
            map_index++;
            registerBehavior( "custom:" + varName, map_index );
            behavior_to_int.put( "custom: " + varName, map_index );
            behavior_to_int.put( "custom " + nc, map_index );
            nc++;
        }

        map_index++;
        registerBehavior( "is_movable", map_index );
        behavior_to_int.put( "movable", map_index );
        behavior_to_int.put( "is movable", map_index );

        // immunogenicity to each cell type 
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            map_index++;
            registerBehavior( "immunogenicity to " + cd.name, map_index );
            behavior_to_int.put( "immunogenicity to cell type " + cd.type, map_index );
        }

        map_index++;
        registerBehavior( "cell attachment rate", map_index );

        map_index++;
        registerBehavior( "cell detachment rate", map_index );

        map_index++;
        registerBehavior( "maximum number of cell attachments", map_index );

        map_index++;
        registerBehavior( "damage rate", map_index );

        /* add new behaviors above this line */

        // resize scales; 
        //        signal_scales.resize( int_to_signal.size() , 1.0 ); 
        signal_scales = VectorUtil.resize( signal_scales, int_to_signal.size(), 1.0 );
        //        System.out.println();
        //        display_signal_dictionary(); 
        //        display_behavior_dictionary(); 
        /*
        // now create empty SR models for each cell definition 
        
        for( int i=0 ; i < cell_definitions_by_index.size() ; i++ )
        { create_SR_model( *cell_definitions_by_index[i] ); }
        */
    }

    static int first_secretion_index;
    static int first_secretion_target_index;
    static int first_uptake_index;
    static int first_export_index;
    static int first_cycle_index;
    static int apoptosis_model_index;
    static int apoptosis_parameter_index;
    static int necrosis_model_index;
    static int necrosis_parameter_index;
    static int migration_speed_index;
    static int migration_bias_index;
    static int persistence_time_index;
    static int first_chemotaxis_index;
    static int cca_index;
    static int elastic_index;
    static int first_affinity_index;
    static int max_adh_distance_index;
    static int ccr_index;
    static int cba_index;
    static int cbr_index;
    static int dead_phago_index;
    static int first_phagocytosis_index;
    static int first_attack_index;
    static int first_fusion_index;
    static int first_transformation_index;
    static int first_custom_ind;
    static int max_custom_ind;
    static int movable_ind;
    static int first_immunogenicity_index;
    static int attachment_rate_ind;
    static int detachment_rate_ind;
    static int max_attachments_ind;
    static int damage_rate_ind;

    static void setSingleBehavior(Cell pCell, int index, double parameter) throws Exception
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();
        int m = microenvironment.number_of_densities();
        int n = CellDefinition.getDefinitionsCount();

        if( index < 0 )
        {
            throw new Exception(
                    "Warning! Tried to set behavior of unknown index " + index + "!" + "         I'll ignore it, but you should fix it!" );
        }

        // substrate-related behaviors 

        // first m entries are secretion 
        first_secretion_index = findBehaviorIndex( microenvironment.density_names[0] + " secretion" ); // 0; 
        if( index >= first_secretion_index && index < first_secretion_index + m )
        {
            pCell.phenotype.secretion.secretionRates[index - first_secretion_index] = parameter;
            return;
        }

        // next m entries are secretion targets
        first_secretion_target_index = findBehaviorIndex( microenvironment.density_names[0] + " secretion target" ); // m; 
        if( index >= first_secretion_target_index && index < first_secretion_target_index + m )
        {
            pCell.phenotype.secretion.saturationDensities[index - first_secretion_target_index] = parameter;
            return;
        }

        // next m entries are uptake rates
        first_uptake_index = findBehaviorIndex( microenvironment.density_names[0] + " uptake" ); // 2*m; 
        if( index >= first_uptake_index && index < first_uptake_index + m )
        {
            pCell.phenotype.secretion.uptakeRates[index - first_uptake_index] = parameter;
            return;
        }

        // next m entries are net export rates 
        first_export_index = findBehaviorIndex( microenvironment.density_names[0] + " export" ); //  3*m; 
        if( index >= first_export_index && index < first_export_index + m )
        {
            pCell.phenotype.secretion.netExportRates[index - first_export_index] = parameter;
            return;
        }

        // cycle entry (exit from phase 0) and exit from up to 5 more phases 
        first_cycle_index = findBehaviorIndex( "exit from cycle phase 0" ); //  4*m; 
        if( index >= first_cycle_index && index < first_cycle_index + 6 )
        {
            int max_cycle_index = pCell.phenotype.cycle.phases.size();
            if( index < first_cycle_index + max_cycle_index )
            {
                pCell.phenotype.cycle.data.setExitRate( index - first_cycle_index, parameter );

                //                pCell.phenotype.cycle.data.exit_rate( index - first_cycle_index ) = parameter;
                return;
            }
            System.out.println( "Warning: Attempted to set a cycle exit rate outside the bounds of the cell's cycle model"
                    + "         Ignoring it, but you should fix this." );
            return;
        }

        // death rates 

        // apoptosis
        apoptosis_model_index = pCell.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );
        apoptosis_parameter_index = findBehaviorIndex( "apoptosis" );
        if( index == apoptosis_parameter_index )
        {
            pCell.phenotype.death.rates.set( apoptosis_model_index, parameter );
            return;
        }

        // necrosis 
        necrosis_model_index = pCell.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
        int necrosis_parameter_index = findBehaviorIndex( "necrosis" );
        if( index == necrosis_parameter_index )
        {
            pCell.phenotype.death.rates.set( necrosis_model_index, parameter );
            return;
        }

        // migration speed
        migration_speed_index = findBehaviorIndex( "migration speed" );
        if( index == migration_speed_index )
        {
            pCell.phenotype.motility.migration_speed = parameter;
            return;
        }

        // migration bias 
        migration_bias_index = findBehaviorIndex( "migration bias" );
        if( index == migration_bias_index )
        {
            pCell.phenotype.motility.migration_bias = parameter;
            return;
        }

        // migration persistence time
        persistence_time_index = findBehaviorIndex( "migration persistence time" );
        if( index == persistence_time_index )
        {
            pCell.phenotype.motility.persistence_time = parameter;
            return;
        }

        // chemotactic sensitivities 
        first_chemotaxis_index = findBehaviorIndex( "chemotactic response to " + microenvironment.density_names[0] );
        if( index >= first_chemotaxis_index && index < first_chemotaxis_index + m )
        {
            pCell.phenotype.motility.chemotactic_sensitivities[index - first_chemotaxis_index] = parameter;
            return;
        }

        // cell-cell adhesion 
        cca_index = findBehaviorIndex( "cell-cell adhesion" );
        if( index == cca_index )
        {
            pCell.phenotype.mechanics.cell_cell_adhesion_strength = parameter;
            return;
        }

        // cell-cell "springs"
        elastic_index = findBehaviorIndex( "cell-cell adhesion elastic constant" );
        if( index == elastic_index )
        {
            pCell.phenotype.mechanics.attachment_elastic_constant = parameter;
            return;
        }

        // cell adhesion affinities 
        first_affinity_index = findBehaviorIndex( "adhesive affinity to " + CellDefinition.getCellDefinition( 0 ).name );
        if( index >= first_affinity_index && index < first_affinity_index + n )
        {
            pCell.phenotype.mechanics.cell_adhesion_affinities[index - first_affinity_index] = parameter;
            return;
        }

        // max relative maximum adhesion distance 
        max_adh_distance_index = findBehaviorIndex( "relative maximum adhesion distance" );
        if( index == max_adh_distance_index )
        {
            pCell.phenotype.mechanics.relative_maximum_adhesion_distance = parameter;
            return;
        }

        // cell-cell repulsion 
        ccr_index = findBehaviorIndex( "cell-cell repulsion" );
        if( index == ccr_index )
        {
            pCell.phenotype.mechanics.cell_cell_repulsion_strength = parameter;
            return;
        }

        // cell-BM adhesion 
        cba_index = findBehaviorIndex( "cell-BM adhesion" );
        if( index == cba_index )
        {
            pCell.phenotype.mechanics.cell_BM_adhesion_strength = parameter;
            return;
        }

        // cell-BM repulsion 
        cbr_index = findBehaviorIndex( "cell-BM repulsion" );
        if( index == cbr_index )
        {
            pCell.phenotype.mechanics.cell_BM_repulsion_strength = parameter;
            return;
        }

        // dead cell phagocytosis
        dead_phago_index = findBehaviorIndex( "phagocytose dead cell" );
        if( index == dead_phago_index )
        {
            pCell.phenotype.cell_interactions.deadPhagocytosisRate = parameter;
            return;
        }

        // phagocytosis of each live cell type 
        first_phagocytosis_index = findBehaviorIndex( "phagocytose " + CellDefinition.getCellDefinition( 0 ).name );
        if( index >= first_phagocytosis_index && index < first_phagocytosis_index + n )
        {
            pCell.phenotype.cell_interactions.livePhagocytosisRates[index - first_phagocytosis_index] = parameter;
            return;
        }

        // attack of each live cell type 
        first_attack_index = findBehaviorIndex( "attack " + CellDefinition.getCellDefinition( 0 ).name );
        if( index >= first_attack_index && index < first_attack_index + n )
        {
            pCell.phenotype.cell_interactions.attackRates[index - first_attack_index] = parameter;
            return;
        }

        // fusion 
        first_fusion_index = findBehaviorIndex( "fuse to " + CellDefinition.getCellDefinition( 0 ).name );
        if( index >= first_fusion_index && index < first_fusion_index + n )
        {
            pCell.phenotype.cell_interactions.fusionRates[index - first_fusion_index] = parameter;
            return;
        }

        // transformation 
        first_transformation_index = findBehaviorIndex( "transform to " + CellDefinition.getCellDefinition( 0 ).name );
        if( index >= first_transformation_index && index < first_transformation_index + n )
        {
            pCell.phenotype.cell_transformations.transformation_rates[index - first_transformation_index] = parameter;
            return;
        }

        // custom behavior
        first_custom_ind = findBehaviorIndex( "custom 0" );
        max_custom_ind = first_custom_ind + pCell.custom_data.variables.size();
        if( first_custom_ind >= 0 && index >= first_custom_ind && index < max_custom_ind )
        {
            pCell.custom_data.variables.get( index - first_custom_ind ).value = parameter;
        }

        // set cell to movable / not movable 
        movable_ind = findBehaviorIndex( "is_movable" );
        if( index == movable_ind )
        {
            pCell.isMovable = parameter > 0.5;
        }

        // immunogenicity to each cell type 
        first_immunogenicity_index = findBehaviorIndex( "immunogenicity to " + CellDefinition.getCellDefinition( 0 ).name );
        if( index >= first_immunogenicity_index && index < first_immunogenicity_index + n )
        {
            pCell.phenotype.cell_interactions.immunogenicities[index - first_immunogenicity_index] = parameter;
            return;
        }

        // set cell attachment rate  
        attachment_rate_ind = findBehaviorIndex( "cell attachment rate" );
        if( index == attachment_rate_ind )
        {
            pCell.phenotype.mechanics.attachment_rate = parameter;
        }

        // set cell detachment rate  
        detachment_rate_ind = findBehaviorIndex( "cell detachment rate" );
        if( index == detachment_rate_ind )
        {
            pCell.phenotype.mechanics.detachment_rate = parameter;
        }

        // maximum number of cell attachments 
        max_attachments_ind = findBehaviorIndex( "maximum number of cell attachments" );
        if( index == max_attachments_ind )
        {
            pCell.phenotype.mechanics.maximum_number_of_attachments = (int)parameter;
        }

        // cell damage rate (for effector attack)
        damage_rate_ind = findBehaviorIndex( "damage rate" );
        if( index == damage_rate_ind )
        {
            pCell.phenotype.cell_interactions.damageRate = parameter;
        }
    }

    public static void setSingleBehavior(Cell cell, String name, double parameter) throws Exception
    {
        int index = findBehaviorIndex( name );
        setSingleBehavior( cell, index, parameter );
    }

    static int findParameterIndex(String responseName)
    {
        if( behavior_to_int.containsKey( responseName ) )
            return behavior_to_int.get( responseName );
        return -1;
    }

    static int findBehaviorIndex(String response_name)
    {
        return findParameterIndex( response_name );
    }
}
