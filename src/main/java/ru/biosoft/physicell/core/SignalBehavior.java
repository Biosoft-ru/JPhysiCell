package ru.biosoft.physicell.core;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;

public class SignalBehavior
{
    public Map<String, Integer> behaviorToInt = new HashMap<>();
    public Map<Integer, String> indexToBehavior = new HashMap<>();
    public Map<Integer, String> indexToSignal = new HashMap<>();
    public Map<String, Integer> signalToIndex = new HashMap<>();
    public double[] signalScales = new double[0];

    public double getSingleSignal(Cell pCell, int index)
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();
        int m = microenvironment.numberDensities();
        int n = pCell.getModel().getDefinitionsCount();

        double out = 0.0;
        if( index < 0 )
        {
            //            std::cout<< "Why would you ask for array[-1]? Why? WHY???? That's it, I quit." << std::endl; 
            return -9e9;
        }

        // first m entries: extracellular concentration 
        int start_substrate_ind = findSignalIndex( microenvironment.densityNames[0] );
        if( start_substrate_ind <= index && index < start_substrate_ind + m )
        {
            out = pCell.nearest_density_vector()[index - start_substrate_ind];
            out /= signalScales[index];
            return out;
        }

        // second m entries: intracellular concentration 
        int start_int_substrate_ind = findSignalIndex( "intracellular " + microenvironment.densityNames[0] );
        if( start_int_substrate_ind <= index && index < start_int_substrate_ind + m )
        {
            out = pCell.phenotype.molecular.internSubstrates[index - start_int_substrate_ind];
            out /= pCell.phenotype.volume.total;
            out /= signalScales[index];
            return out;
        }

        // next m entries: gradients 
        int start_substrate_grad_ind = findSignalIndex( microenvironment.densityNames[0] + " gradient" );
        if( start_substrate_grad_ind <= index && index < start_substrate_grad_ind + m )
        {
            out = VectorUtil.norm( pCell.nearest_gradient( index - start_substrate_grad_ind ) );
            out /= signalScales[index];
            return out;
        }

        // mechanical pressure 
        int pressure_ind = findSignalIndex( "pressure" );
        if( index == pressure_ind )
        {
            out = pCell.state.simplePressure;
            out /= signalScales[index];
            return out;
        }

        // cell volume  
        int volume_ind = findSignalIndex( "volume" );
        if( index == volume_ind )
        {
            out = pCell.phenotype.volume.total;
            out /= signalScales[index];
            return out;
        }

        // physical contact with cells (of each type) 
        // individual contact signals are a bit costly 
        int contact_ind = findSignalIndex( "contact with " + pCell.getModel().getCellDefinition( 0 ).name );
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
                out /= signalScales[index];
                return out;
            }

            int live_contact_ind = findSignalIndex( "contact with live cell" );
            if( index == live_contact_ind )
            {
                out = live_cells;
                out /= signalScales[index];
                return out;
            }

            int dead_contact_ind = findSignalIndex( "contact with dead cell" );
            // index == dead_contact_ind TODO: check
            out = dead_cells;
            out /= signalScales[index];
            return out;
        }

        // contact with BM 
        int BM_contact_ind = findSignalIndex( "contact with basement membrane" );
        if( index == BM_contact_ind )
        {
            //            out = (double) pCell.state.contact_with_basement_membrane; 
            out = pCell.state.contactWithBasementMembrane ? 1 : 0;
            out /= signalScales[index];
            return out;
        }

        // damage
        int damage_ind = findSignalIndex( "damage" );
        if( index == damage_ind )
        {
            out = pCell.state.damage;
            out /= signalScales[index];
            return out;
        }

        // live / dead state 
        int dead_ind = findSignalIndex( "dead" );
        if( index == dead_ind )
        {
            out = pCell.phenotype.death.dead ? 1 : 0;
            out /= signalScales[index];
            return out;
        }

        // integrated total attack time 
        int tot_attack_ind = findSignalIndex( "total attack time" );
        if( index == tot_attack_ind )
        {
            out = pCell.state.totalAttackTime;
            out /= signalScales[index];
            return out;
        }

        // time 
        int time_ind = findSignalIndex( "time" );
        if( index == time_ind )
        {
            out = pCell.getMicroenvironment().time;
            out /= signalScales[index]; //TODO: later
            return out;
        }

        // custom signals 
        int first_custom_ind = findSignalIndex( "custom 0" );
        int max_custom_ind = first_custom_ind + pCell.customData.variables.size();
        if( first_custom_ind > -1 && index >= first_custom_ind && index < max_custom_ind )
        {
            out = pCell.customData.variables.get( index - first_custom_ind ).value;
            out /= signalScales[index];
            return out;
        }

        int apoptotic_ind = findSignalIndex( "apoptotic" );
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

        int necrotic_ind = findSignalIndex( "necrotic" );
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

    public static Set<String> getSignals(String[] densities, String[] cellTypes, String[] custom)
    {
        Set<String> signals = new HashSet<>();
        signals.add( "pressure" );
        signals.add( "volume" );
        signals.add( "contact with live cell" );
        signals.add( "contact with dead cell" );
        signals.add( "contact with basement membrane" );
        signals.add( "damage" );
        signals.add( "is dead" );
        signals.add( "total attack time" );
        signals.add( "time" );
        signals.add( "apoptotic" );
        signals.add( "nectrotic" );
        for( String s : densities )
        {
            signals.add( s );
            signals.add( "intracellular " + s );
            signals.add( "gradient of " + s );
        }
        Stream.of( cellTypes ).forEach( s -> signals.add( "contact with " + s ) );
        Stream.of( custom ).forEach( s -> signals.add( "custom: " + s ) );
        return signals;
    }

    public static Set<String> getBehaviors(String[] densities, String[] cellTypes, String[] custom)
    {
        Set<String> behaviors = new HashSet<String>();
        behaviors.add( "cycle entry" );
        behaviors.add( "exit from cycle phase " );//number of phases
        behaviors.add( "apoptosis" );
        behaviors.add( "necrosis" );
        behaviors.add( "apoptosis" );
        behaviors.add( "migration speed" );
        behaviors.add( "migration bias" );
        behaviors.add( "migration persistence time" );
        behaviors.add( "cell-cell adhesion" );
        behaviors.add( "cell-cell adhesion elastic constant" );
        behaviors.add( "relative maximum adhesion distance" );
        behaviors.add( "cell-cell repulsion" );
        behaviors.add( "cell-membrane adhesion" );
        behaviors.add( "cell-membrane repulsion" );
        behaviors.add( "phagocytosis of dead cells" );
        behaviors.add( "cell attachment rate" );
        behaviors.add( "cell detachment rate" );
        behaviors.add( "maximum number of cell attachments" );
        behaviors.add( "damage rate" );
        behaviors.add( "is movable" );

        for( String cellType : cellTypes )
        {
            behaviors.add( "adhesive affinity to " + cellType );
            behaviors.add( "phagocytosis of " + cellType );
            behaviors.add( "attack " + cellType );
            behaviors.add( "fuse to " + cellType );
            behaviors.add( "transform to " + cellType );
            behaviors.add( "immunogenicity to " + cellType );
        }

        for( String d : densities )
        {
            behaviors.add( d + " secretion" );
            behaviors.add( d + " secretion target" );
            behaviors.add( d + " secretion saturation density" );
            behaviors.add( d + " uptake" );
            behaviors.add( d + " export" );
            behaviors.add( "chemotactic response to " + d );
        }
        Stream.of( custom ).forEach( s -> behaviors.add( "custom: " + s ) );

        return behaviors;
    }

    int findSignalIndex(String name)
    {
        if( signalToIndex.containsKey( name ) )
            return signalToIndex.get( name );
        System.out.println( "having trouble finding " + name );
        return -1;
    }

    public double getSingleSignal(Cell pCell, String name)
    {
        return getSingleSignal( pCell, findSignalIndex( name ) );
    }

    private void registerSignal(String name, int index)
    {
        signalToIndex.put( name, index );// ] = m+i;
        indexToSignal.put( index, name );//] = name; 
    }

    private void registerBehavior(String name, int index)
    {
        behaviorToInt.put( name, index );// ] = m+i;
        indexToBehavior.put( index, name );//] = name; 
    }

    public void setupDictionaries(Model model)
    {
        Microenvironment microenvironment = model.getMicroenvironment();
        int m = microenvironment.numberDensities();

        signalToIndex.clear();
        indexToSignal.clear();

        for( int i = 0; i < m; i++ )
        {
            registerSignal( microenvironment.densityNames[i], i );
        }

        // internalized substrates 
        int mapIndex = m;
        for( int i = 0; i < m; i++ )
        {
            registerSignal( "intracellular " + microenvironment.densityNames[i], m + i );
            signalToIndex.put( "internalized " + microenvironment.densityNames[i], m + i );//[ name ] = m+i; 
        }

        // substrate gradients 
        mapIndex = 2 * m;
        for( int i = 0; i < m; i++ )
        {
            registerSignal( microenvironment.densityNames[i] + " gradient", mapIndex );
            signalToIndex.put( "grad(" + microenvironment.densityNames[i] + ")", mapIndex );//[ name ] = map_index;
            signalToIndex.put( "gradient of " + microenvironment.densityNames[i], mapIndex );// ] = map_index;
            mapIndex++;
        }

        // mechanical pressure 
        registerSignal( "pressure", mapIndex );

        // total volume 
        mapIndex++;
        registerSignal( "volume", mapIndex );

        // contact with each cell type 
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            mapIndex++;
            signalToIndex.put( "contact with " + cd.name, mapIndex );
            registerSignal( "contact with cell type " + cd.type, mapIndex );
        }

        // contact with (any) live cell 
        mapIndex++;
        registerSignal( "contact with live cell", mapIndex );
        signalToIndex.put( "contact with live cells", mapIndex );

        // contact with dead cell 
        mapIndex++;
        registerSignal( "contact with dead cell", mapIndex );
        signalToIndex.put( "contact with dead cells", mapIndex );

        // contact with basement membrane 
        mapIndex++;
        registerSignal( "contact with basement membrane", mapIndex );
        signalToIndex.put( "contact with BM", mapIndex );

        // damage state 
        mapIndex++;
        registerSignal( "damage", mapIndex );

        // live / dead state 
        mapIndex++;
        registerSignal( "dead", mapIndex );
        signalToIndex.put( "is dead", mapIndex );

        // total attack time 
        mapIndex++;
        registerSignal( "total attack time", mapIndex );

        // current time
        mapIndex++;
        registerSignal( "time", mapIndex );
        signalToIndex.put( "current time", mapIndex );
        signalToIndex.put( "global time", mapIndex );

        int nc = 0;
        //        Set<String> customVars = new HashSet<>();
        //        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        //        {
        //            for( Variable var : cd.custom_data.variables )
        //                customVars.add( var.name );
        //        }
        CellDefinition def = model.getCellDefinition( 0 );
        // custom signals
        for( Variable var : def.custom_data.variables )
        {
            String varName = var.name;
            mapIndex++;
            registerSignal( "custom:" + varName, mapIndex );
            signalToIndex.put( "custom: " + varName, mapIndex );
            signalToIndex.put( "custom " + nc, mapIndex );
            nc++;
        }

        mapIndex++;
        registerSignal( "apoptotic", mapIndex );
        signalToIndex.put( "is_apoptotic", mapIndex );

        mapIndex++;
        registerSignal( "necrotic", mapIndex );
        signalToIndex.put( "is_necrotic", mapIndex );

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
        behaviorToInt.clear();
        indexToBehavior.clear();

        // construct behaviors 
        String name;
        for( int i = 0; i < m; i++ )
        {
            mapIndex = i;
            name = microenvironment.densityNames[i];

            // secretion rate 
            registerBehavior( name + " " + "secretion", mapIndex );

            // secretion target 
            mapIndex = m + i;
            registerBehavior( name + " " + "secretion target", mapIndex );

            // synonym 
            behaviorToInt.put( name + " " + "secretion saturation density", mapIndex );

            // uptake rate 
            mapIndex = 2 * m + i;
            registerBehavior( name + " " + "uptake", mapIndex );

            // net export rate 
            mapIndex = 3 * m + i;
            registerBehavior( name + " " + "export", mapIndex );
        }

        mapIndex = 4 * m;
        registerBehavior( "cycle entry", mapIndex );

        // synonym 
        behaviorToInt.put( "exit from cycle phase 0", mapIndex );

        // other cyle phases 
        for( int i = 1; i < 6; i++ )
        {
            mapIndex++;
            registerBehavior( "exit from cycle phase " + i, mapIndex );
        }

        mapIndex++;
        registerBehavior( "apoptosis", mapIndex );

        mapIndex++;
        registerBehavior( "necrosis", mapIndex );

        mapIndex++;
        registerBehavior( "migration speed", mapIndex );

        mapIndex++;
        registerBehavior( "migration bias", mapIndex );

        mapIndex++;
        registerBehavior( "migration persistence time", mapIndex );

        // chemotactic sensitivities 
        for( int i = 0; i < m; i++ )
        {
            mapIndex++;
            registerBehavior( "chemotactic response to " + microenvironment.densityNames[i], mapIndex );
            behaviorToInt.put( "chemotactic sensitivity to " + microenvironment.densityNames[i], mapIndex );
        }

        // cell-cell adhesion 
        mapIndex++;
        registerBehavior( "cell-cell adhesion", mapIndex );

        mapIndex++;
        registerBehavior( "cell-cell adhesion elastic constant", mapIndex );

        // cell adhesion affinities 
        // cell-type specific adhesion 
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            mapIndex++;
            registerBehavior( "adhesive affinity to " + cd.name, mapIndex );
            behaviorToInt.put( "adhesive affinity to cell type " + cd.type, mapIndex );
        }

        // max adhesion distance 
        mapIndex++;
        registerBehavior( "relative maximum adhesion distance", mapIndex );

        // cell-cell repulsion 
        mapIndex++;
        registerBehavior( "cell-cell repulsion", mapIndex );

        // cell-BM adhesion 
        mapIndex++;
        registerBehavior( "cell-BM adhesion", mapIndex );
        behaviorToInt.put( "cell-membrane adhesion", mapIndex );

        // cell-BM repulsion 
        mapIndex++;
        registerBehavior( "cell-BM repulsion", mapIndex );
        behaviorToInt.put( "cell-membrane repulsion", mapIndex );

        // phagocytosis of dead cell
        mapIndex++;
        registerBehavior( "phagocytose dead cell", mapIndex );
        behaviorToInt.put( "phagocytosis of dead cell", mapIndex );
        behaviorToInt.put( "phagocytosis of dead cells", mapIndex );

        // phagocytosis of each live cell type 
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            mapIndex++;
            registerBehavior( "phagocytose " + cd.name, mapIndex );
            behaviorToInt.put( "phagocytose cell type " + cd.type, mapIndex );
            behaviorToInt.put( "phagocytosis of " + cd.type, mapIndex );
        }

        // attack of each live cell type 
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            mapIndex++;
            registerBehavior( "attack " + cd.name, mapIndex );
            behaviorToInt.put( "attack cell type " + cd.type, mapIndex );
        }

        // fusion 
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            mapIndex++;
            registerBehavior( "fuse to " + cd.name, mapIndex );
            behaviorToInt.put( "fuse to cell type " + cd.type, mapIndex );
        }

        // transformation 
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            mapIndex++;
            registerBehavior( "transform to " + cd.name, mapIndex );
            behaviorToInt.put( "transform to cell type " + cd.type, mapIndex );
        }

        // custom behaviors
        nc = 0;
        for( Variable var : def.custom_data.variables )
        {
            String varName = var.name;
            mapIndex++;
            registerBehavior( "custom:" + varName, mapIndex );
            behaviorToInt.put( "custom: " + varName, mapIndex );
            behaviorToInt.put( "custom " + nc, mapIndex );
            nc++;
        }

        mapIndex++;
        registerBehavior( "is_movable", mapIndex );
        behaviorToInt.put( "movable", mapIndex );
        behaviorToInt.put( "is movable", mapIndex );

        // immunogenicity to each cell type 
        for( CellDefinition cd : model.getCellDefinitions() )
        {
            mapIndex++;
            registerBehavior( "immunogenicity to " + cd.name, mapIndex );
            behaviorToInt.put( "immunogenicity to cell type " + cd.type, mapIndex );
        }

        mapIndex++;
        registerBehavior( "cell attachment rate", mapIndex );

        mapIndex++;
        registerBehavior( "cell detachment rate", mapIndex );

        mapIndex++;
        registerBehavior( "maximum number of cell attachments", mapIndex );

        mapIndex++;
        registerBehavior( "damage rate", mapIndex );

        /* add new behaviors above this line */

        // resize scales; 
        //        signal_scales.resize( int_to_signal.size() , 1.0 ); 
        signalScales = VectorUtil.resize( signalScales, indexToSignal.size(), 1.0 );
        //        System.out.println();
        //        display_signal_dictionary(); 
        //        display_behavior_dictionary(); 
        /*
        // now create empty SR models for each cell definition 
        
        for( int i=0 ; i < cell_definitions_by_index.size() ; i++ )
        { create_SR_model( *cell_definitions_by_index[i] ); }
        */
    }

    //    static int first_secretion_index;
    //    static int first_secretion_target_index;
    //    static int first_uptake_index;
    //    static int first_export_index;
    //    static int first_cycle_index;
    //    static int apoptosis_model_index;
    //    static int apoptosis_parameter_index;
    //    static int necrosis_model_index;
    //    static int necrosis_parameter_index;
    //    static int migration_speed_index;
    //    static int migration_bias_index;
    //    static int persistence_time_index;
    //    static int first_chemotaxis_index;
    //    static int cca_index;
    //    static int elastic_index;
    //    static int first_affinity_index;
    //    static int max_adh_distance_index;
    //    static int ccr_index;
    //    static int cba_index;
    //    static int cbr_index;
    //    static int dead_phago_index;
    //    static int first_phagocytosis_index;
    //    static int first_attack_index;
    //    static int first_fusion_index;
    //    static int first_transformation_index;
    //    static int first_custom_ind;
    //    static int max_custom_ind;
    //    static int movable_ind;
    //    static int first_immunogenicity_index;
    //    static int attachment_rate_ind;
    //    static int detachment_rate_ind;
    //    static int max_attachments_ind;
    //    static int damage_rate_ind;

    void setSingleBehavior(Cell pCell, int index, double parameter) throws Exception
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();
        Model model = pCell.getModel();
        int m = microenvironment.numberDensities();
        int n = model.getDefinitionsCount();

        if( index < 0 )
        {
            throw new Exception(
                    "Warning! Tried to set behavior of unknown index " + index + "!" + "         I'll ignore it, but you should fix it!" );
        }

        // substrate-related behaviors 

        // first m entries are secretion 
        int first_secretion_index = findBehaviorIndex( microenvironment.densityNames[0] + " secretion" ); // 0; 
        if( index >= first_secretion_index && index < first_secretion_index + m )
        {
            pCell.phenotype.secretion.secretionRates[index - first_secretion_index] = parameter;
            return;
        }

        // next m entries are secretion targets
        int first_secretion_target_index = findBehaviorIndex( microenvironment.densityNames[0] + " secretion target" ); // m; 
        if( index >= first_secretion_target_index && index < first_secretion_target_index + m )
        {
            pCell.phenotype.secretion.saturationDensities[index - first_secretion_target_index] = parameter;
            return;
        }

        // next m entries are uptake rates
        int first_uptake_index = findBehaviorIndex( microenvironment.densityNames[0] + " uptake" ); // 2*m; 
        if( index >= first_uptake_index && index < first_uptake_index + m )
        {
            pCell.phenotype.secretion.uptakeRates[index - first_uptake_index] = parameter;
            return;
        }

        // next m entries are net export rates 
        int first_export_index = findBehaviorIndex( microenvironment.densityNames[0] + " export" ); //  3*m; 
        if( index >= first_export_index && index < first_export_index + m )
        {
            pCell.phenotype.secretion.netExportRates[index - first_export_index] = parameter;
            return;
        }

        // cycle entry (exit from phase 0) and exit from up to 5 more phases 
        int first_cycle_index = findBehaviorIndex( "exit from cycle phase 0" ); //  4*m; 
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
        int apoptosis_model_index = pCell.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );
        int apoptosis_parameter_index = findBehaviorIndex( "apoptosis" );
        if( index == apoptosis_parameter_index )
        {
            pCell.phenotype.death.rates.set( apoptosis_model_index, parameter );
            return;
        }

        // necrosis 
        int necrosis_model_index = pCell.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
        int necrosis_parameter_index = findBehaviorIndex( "necrosis" );
        if( index == necrosis_parameter_index )
        {
            pCell.phenotype.death.rates.set( necrosis_model_index, parameter );
            return;
        }

        // migration speed
        int migration_speed_index = findBehaviorIndex( "migration speed" );
        if( index == migration_speed_index )
        {
            pCell.phenotype.motility.migrationSpeed = parameter;
            return;
        }

        // migration bias 
        int migration_bias_index = findBehaviorIndex( "migration bias" );
        if( index == migration_bias_index )
        {
            pCell.phenotype.motility.migrationBias = parameter;
            return;
        }

        // migration persistence time
        int persistence_time_index = findBehaviorIndex( "migration persistence time" );
        if( index == persistence_time_index )
        {
            pCell.phenotype.motility.persistenceTime = parameter;
            return;
        }

        // chemotactic sensitivities 
        int first_chemotaxis_index = findBehaviorIndex( "chemotactic response to " + microenvironment.densityNames[0] );
        if( index >= first_chemotaxis_index && index < first_chemotaxis_index + m )
        {
            pCell.phenotype.motility.chemotacticSensitivities[index - first_chemotaxis_index] = parameter;
            return;
        }

        // cell-cell adhesion 
        int cca_index = findBehaviorIndex( "cell-cell adhesion" );
        if( index == cca_index )
        {
            pCell.phenotype.mechanics.cellCellAdhesionStrength = parameter;
            return;
        }

        // cell-cell "springs"
        int elastic_index = findBehaviorIndex( "cell-cell adhesion elastic constant" );
        if( index == elastic_index )
        {
            pCell.phenotype.mechanics.attachmentElasticConstant = parameter;
            return;
        }

        // cell adhesion affinities 
        int first_affinity_index = findBehaviorIndex( "adhesive affinity to " + pCell.getModel().getCellDefinition( 0 ).name );
        if( index >= first_affinity_index && index < first_affinity_index + n )
        {
            pCell.phenotype.mechanics.cellAdhesionAffinities[index - first_affinity_index] = parameter;
            return;
        }

        // max relative maximum adhesion distance 
        int max_adh_distance_index = findBehaviorIndex( "relative maximum adhesion distance" );
        if( index == max_adh_distance_index )
        {
            pCell.phenotype.mechanics.relMaxAdhesionDistance = parameter;
            return;
        }

        // cell-cell repulsion 
        int ccr_index = findBehaviorIndex( "cell-cell repulsion" );
        if( index == ccr_index )
        {
            pCell.phenotype.mechanics.cellCellRepulsionStrength = parameter;
            return;
        }

        // cell-BM adhesion 
        int cba_index = findBehaviorIndex( "cell-BM adhesion" );
        if( index == cba_index )
        {
            pCell.phenotype.mechanics.cellBMAdhesionStrength = parameter;
            return;
        }

        // cell-BM repulsion 
        int cbr_index = findBehaviorIndex( "cell-BM repulsion" );
        if( index == cbr_index )
        {
            pCell.phenotype.mechanics.cellBMRepulsionStrength = parameter;
            return;
        }

        // dead cell phagocytosis
        int dead_phago_index = findBehaviorIndex( "phagocytose dead cell" );
        if( index == dead_phago_index )
        {
            pCell.phenotype.cellInteractions.deadPhagocytosisRate = parameter;
            return;
        }

        // phagocytosis of each live cell type 
        int first_phagocytosis_index = findBehaviorIndex( "phagocytose " + model.getCellDefinition( 0 ).name );
        if( index >= first_phagocytosis_index && index < first_phagocytosis_index + n )
        {
            pCell.phenotype.cellInteractions.livePhagocytosisRates[index - first_phagocytosis_index] = parameter;
            return;
        }

        // attack of each live cell type 
        int first_attack_index = findBehaviorIndex( "attack " + model.getCellDefinition( 0 ).name );
        if( index >= first_attack_index && index < first_attack_index + n )
        {
            pCell.phenotype.cellInteractions.attackRates[index - first_attack_index] = parameter;
            return;
        }

        // fusion 
        int first_fusion_index = findBehaviorIndex( "fuse to " + model.getCellDefinition( 0 ).name );
        if( index >= first_fusion_index && index < first_fusion_index + n )
        {
            pCell.phenotype.cellInteractions.fusionRates[index - first_fusion_index] = parameter;
            return;
        }

        // transformation 
        int first_transformation_index = findBehaviorIndex( "transform to " + model.getCellDefinition( 0 ).name );
        if( index >= first_transformation_index && index < first_transformation_index + n )
        {
            pCell.phenotype.cellTransformations.transformationRates[index - first_transformation_index] = parameter;
            return;
        }

        // custom behavior
        int first_custom_ind = findBehaviorIndex( "custom 0" );
        int max_custom_ind = first_custom_ind + pCell.customData.variables.size();
        if( first_custom_ind >= 0 && index >= first_custom_ind && index < max_custom_ind )
        {
            pCell.customData.variables.get( index - first_custom_ind ).value = parameter;
        }

        // set cell to movable / not movable 
        int movable_ind = findBehaviorIndex( "is_movable" );
        if( index == movable_ind )
        {
            pCell.isMovable = parameter > 0.5;
        }

        // immunogenicity to each cell type 
        int first_immunogenicity_index = findBehaviorIndex( "immunogenicity to " + model.getCellDefinition( 0 ).name );
        if( index >= first_immunogenicity_index && index < first_immunogenicity_index + n )
        {
            pCell.phenotype.cellInteractions.immunogenicities[index - first_immunogenicity_index] = parameter;
            return;
        }

        // set cell attachment rate  
        int attachment_rate_ind = findBehaviorIndex( "cell attachment rate" );
        if( index == attachment_rate_ind )
        {
            pCell.phenotype.mechanics.attachmentRate = parameter;
        }

        // set cell detachment rate  
        int detachment_rate_ind = findBehaviorIndex( "cell detachment rate" );
        if( index == detachment_rate_ind )
        {
            pCell.phenotype.mechanics.detachmentRate = parameter;
        }

        // maximum number of cell attachments 
        int max_attachments_ind = findBehaviorIndex( "maximum number of cell attachments" );
        if( index == max_attachments_ind )
        {
            pCell.phenotype.mechanics.maxAttachments = (int)parameter;
        }

        // cell damage rate (for effector attack)
        int damage_rate_ind = findBehaviorIndex( "damage rate" );
        if( index == damage_rate_ind )
        {
            pCell.phenotype.cellInteractions.damageRate = parameter;
        }
    }

    public void setSingleBehavior(Cell cell, String name, double parameter) throws Exception
    {
        int index = findBehaviorIndex( name );
        setSingleBehavior( cell, index, parameter );
    }

    public int findBehaviorIndex(String responseName)
    {
        Integer result = behaviorToInt.get( responseName );
        if( result == null )
            return -1;
        return result;
    }

    public double getSinglBehavior(Cell pCell, String name) throws Exception
    {
        return getSingleBehavior( pCell, findBehaviorIndex( name ) );
    }

    public double getSingleBehavior(Cell pCell, int index) throws Exception
    {
        Model model = pCell.getModel();
        Microenvironment microenvironment = pCell.getMicroenvironment();
        int m = microenvironment.numberDensities();
        //        static int m = microenvironment.number_of_densities(); 
        int n = pCell.getModel().getDefinitionsCount();//.size(); 

        if( index < 0 )
        {
            throw new Exception(
                    "Warning: attempted to get behavior with unknown index " + "         I'm ignoring it, but you should fix it!" );
        }

        // substrate-related behaviors 

        // first m entries are secretion 
        int first_secretion_index = findBehaviorIndex( microenvironment.densityNames[0] + " secretion" ); // 0; 
        if( index >= first_secretion_index && index < first_secretion_index + m )
        {
            return pCell.phenotype.secretion.secretionRates[index - first_secretion_index];
        }

        // next m entries are secretion targets
        int first_secretion_target_index = findBehaviorIndex( microenvironment.densityNames[0] + " secretion target" ); // m; 
        if( index >= first_secretion_target_index && index < first_secretion_target_index + m )
        {
            return pCell.phenotype.secretion.saturationDensities[index - first_secretion_target_index];
        }

        // next m entries are uptake rates
        int first_uptake_index = findBehaviorIndex( microenvironment.densityNames[0] + " uptake" ); // 2*m; 
        if( index >= first_uptake_index && index < first_uptake_index + m )
        {
            return pCell.phenotype.secretion.uptakeRates[index - first_uptake_index];
        }

        // next m entries are net export rates 
        int first_export_index = findBehaviorIndex( microenvironment.densityNames[0] + " export" ); //  3*m; 
        if( index >= first_export_index && index < first_export_index + m )
        {
            return pCell.phenotype.secretion.netExportRates[index - first_export_index];
        }

        // cycle entry (exit from phase 0) and exit from up to 5 more phases 
        int first_cycle_index = findBehaviorIndex( "exit from cycle phase 0" ); //  4*m; 
        int max_cycle_index = pCell.phenotype.cycle.phases.size();
        if( max_cycle_index > 6 )
        {
            max_cycle_index = 6;
            System.out.println( "Warning: Standardized behaviors only support exit rate from the first 6 phases of a cell cycle!"
                    + "         Ignoring any later phase exit rates." );
        }
        if( index >= first_cycle_index && index < first_cycle_index + 6 )
        {
            int ind = index - first_cycle_index;
            if( ind < max_cycle_index )
            {
                return pCell.phenotype.cycle.data.getExitRate( ind );
            }//.exit_rate( ind ); }
            return 0.0;
        }

        int apoptosis_index = pCell.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );
        int necrosis_index = pCell.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );

        int apop_param_index = findBehaviorIndex( "apoptosis" );
        int necr_param_index = findBehaviorIndex( "necrosis" );

        // apoptosis 
        if( index == apop_param_index )
        {
            return pCell.phenotype.death.rates.get( apoptosis_index );
        }

        // necrosis 
        if( index == necr_param_index )
        {
            return pCell.phenotype.death.rates.get( necrosis_index );
        }

        // migration speed
        int migr_spd_index = findBehaviorIndex( "migration speed" );
        if( index == migr_spd_index )
        {
            return pCell.phenotype.motility.migrationSpeed;
        }

        // migration bias 
        int migr_bias_index = findBehaviorIndex( "migration bias" );
        if( index == migr_bias_index )
        {
            return pCell.phenotype.motility.migrationBias;
        }

        // migration persistence time
        int migr_pt_index = findBehaviorIndex( "migration persistence time" );
        if( index == migr_pt_index )
        {
            return pCell.phenotype.motility.persistenceTime;
        }

        // chemotactic sensitivities 
        int first_chemotaxis_index = findBehaviorIndex( "chemotactic response to " + microenvironment.densityNames[0] );
        if( index >= first_chemotaxis_index && index < first_chemotaxis_index + m )
        {
            return pCell.phenotype.motility.chemotacticSensitivities[index - first_chemotaxis_index];
        }

        // cell-cell adhesion 
        int cca_index = findBehaviorIndex( "cell-cell adhesion" );
        if( index == cca_index )
        {
            return pCell.phenotype.mechanics.cellCellAdhesionStrength;
        }

        // cell-cell "springs"
        int cca_spring_index = findBehaviorIndex( "cell-cell adhesion elastic constant" );
        if( index == cca_spring_index )
        {
            return pCell.phenotype.mechanics.attachmentElasticConstant;
        }

        // cell adhesion affinities 
        int first_affinity_index = findBehaviorIndex( "adhesive affinity to " + pCell.getModel().getCellDefinition( 0 ).name );//cell_definitions_by_type[0].name ); 
        if( index >= first_affinity_index && index < first_affinity_index + n )
        {
            return pCell.phenotype.mechanics.cellAdhesionAffinities[index - first_affinity_index];
        }

        // max relative maximum adhesion distance 
        int max_adh_index = findBehaviorIndex( "relative maximum adhesion distance" );
        if( index == max_adh_index )
        {
            return pCell.phenotype.mechanics.relMaxAdhesionDistance;
        }

        // cell-cell repulsion 
        int ccr_index = findBehaviorIndex( "cell-cell repulsion" );
        if( index == ccr_index )
        {
            return pCell.phenotype.mechanics.cellCellRepulsionStrength;
        }

        // cell-BM adhesion 
        int cba_index = findBehaviorIndex( "cell-BM adhesion" );
        if( index == cba_index )
        {
            return pCell.phenotype.mechanics.cellBMAdhesionStrength;
        }

        // cell-BM repulsion 
        int cbr_index = findBehaviorIndex( "cell-BM repulsion" );
        if( index == cbr_index )
        {
            return pCell.phenotype.mechanics.cellBMRepulsionStrength;
        }

        // dead cell phagocytosis
        int dead_phag_index = findBehaviorIndex( "phagocytose dead cell" );
        if( index == dead_phag_index )
        {
            return pCell.phenotype.cellInteractions.deadPhagocytosisRate;
        }

        // phagocytosis of each live cell type 
        int first_phagocytosis_index = findBehaviorIndex( "phagocytose " + model.getCellDefinition( 0 ).name );//cell_definitions_by_type[0].name ); 
        if( index >= first_phagocytosis_index && index < first_phagocytosis_index + n )
        {
            return pCell.phenotype.cellInteractions.livePhagocytosisRates[index - first_phagocytosis_index];
        }

        // attack of each live cell type 
        int first_attack_index = findBehaviorIndex( "attack " + model.getCellDefinition( 0 ).name );//cell_definitions_by_type[0].name ); 
        if( index >= first_attack_index && index < first_attack_index + n )
        {
            return pCell.phenotype.cellInteractions.attackRates[index - first_attack_index];
        }

        // fusion 
        int first_fusion_index = findBehaviorIndex( "fuse to " + model.getCellDefinition( 0 ).name );//cell_definitions_by_type[0].name ); 
        if( index >= first_fusion_index && index < first_fusion_index + n )
        {
            return pCell.phenotype.cellInteractions.fusionRates[index - first_fusion_index];
        }

        // transformation 
        int first_transformation_index = findBehaviorIndex( "transform to " + model.getCellDefinition( 0 ).name );//cell_definitions_by_type[0].name ); 
        if( index >= first_transformation_index && index < first_transformation_index + n )
        {
            return pCell.phenotype.cellTransformations.transformationRates[index - first_transformation_index];
        }

        // custom behavior
        int first_custom_ind = findBehaviorIndex( "custom 0" );
        int max_custom_ind = first_custom_ind + pCell.customData.variables.size();
        if( first_custom_ind >= 0 && index >= first_custom_ind && index < max_custom_ind )
        {
            return pCell.customData.variables.get( index - first_custom_ind ).value;
        }

        // is the cell movable / not movable 
        int movable_ind = findBehaviorIndex( "is_movable" );
        if( index == movable_ind )
        {
            if( pCell.isMovable == true )
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        }

        // vector of immunogenicity behaviors 
        int start_immunogenicity_ind = findBehaviorIndex( "immunogenicity to " + model.getCellDefinition( 0 ).name );// Dcell_definitions_by_type[0].name ); 
        int max_immunogenicity_ind = start_immunogenicity_ind + n;
        if( start_immunogenicity_ind > -1 && index >= start_immunogenicity_ind && index < max_immunogenicity_ind )
        {
            return pCell.phenotype.cellInteractions.immunogenicities[index - start_immunogenicity_ind];
        }


        // set cell attachment rate  
        int attachment_rate_ind = findBehaviorIndex( "cell attachment rate" );
        if( index == attachment_rate_ind )
            return pCell.phenotype.mechanics.attachmentRate;

        // set cell detachment rate  
        int detachment_rate_ind = findBehaviorIndex( "cell detachment rate" );
        if( index == detachment_rate_ind )
        {
            return pCell.phenotype.mechanics.detachmentRate;
        }

        // maximum number of cell attachments 
        int max_attachments_ind = findBehaviorIndex( "maximum number of cell attachments" );
        if( index == max_attachments_ind )
        {
            return pCell.phenotype.mechanics.maxAttachments;
        }

        // get damage rate 
        int damage_rate_ind = findBehaviorIndex( "damage rate" );
        if( index == damage_rate_ind )
        {
            return pCell.phenotype.cellInteractions.damageRate;
        }

        return -1;
    }

    public double[] getBaseBehaviors(Cell pCell)
    {
        CellDefinition pCD = pCell.getModel().getCellDefinition( pCell.typeName );
        Model model = pCell.getModel();
        Microenvironment microenvironment = pCell.getMicroenvironment();

        double[] parameters = new double[indexToBehavior.size()];//( int_to_behavior.size() , 0.0 ); 

        // substrate-related behaviors 

        // first m entries are secretion 
        int first_secretion_index = findBehaviorIndex( microenvironment.densityNames[0] + " secretion" ); // 0; 
        System.arraycopy( pCD.phenotype.secretion.secretionRates, 0, parameters, first_secretion_index,
                pCD.phenotype.secretion.secretionRates.length );
        //        std::copy(  pCD.phenotype.secretion.secretionRates.begin(), 
        //                    pCD.phenotype.secretion.secretionRates.end(), 
        //                    parameters.begin()+first_secretion_index ); 


        // next m entries are secretion targets
        int first_secretion_target_index = findBehaviorIndex( microenvironment.densityNames[0] + " secretion target" ); // m; 
        System.arraycopy( pCD.phenotype.secretion.saturationDensities, 0, parameters, first_secretion_target_index,
                pCD.phenotype.secretion.saturationDensities.length );
        //        std::copy(  pCD.phenotype.secretion.saturation_densities.begin(), 
        //                    pCD.phenotype.secretion.saturation_densities.end(), 
        //                    parameters.begin()+first_secretion_target_index ); 

        // next m entries are uptake rates
        int first_uptake_index = findBehaviorIndex( microenvironment.densityNames[0] + " uptake" ); // 2*m; 
        System.arraycopy( pCD.phenotype.secretion.uptakeRates, 0, parameters, first_uptake_index,
                pCD.phenotype.secretion.uptakeRates.length );
        //        std::copy(  pCD.phenotype.secretion.uptake_rates.begin(), 
        //                    pCD.phenotype.secretion.uptake_rates.end(), 
        //                    parameters.begin()+first_uptake_index ); 

        // next m entries are net export rates 
        int first_export_index = findBehaviorIndex( microenvironment.densityNames[0] + " export" ); //  3*m; 
        System.arraycopy( pCD.phenotype.secretion.netExportRates, 0, parameters, first_export_index,
                pCD.phenotype.secretion.netExportRates.length );
        //        std::copy(  pCD.phenotype.secretion.net_export_rates.begin(), 
        //                    pCD.phenotype.secretion.net_export_rates.end(), 
        //                    parameters.begin()+first_export_index ); 

        // cycle entry (exit from phase 0) and exit from up to 5 more phases 
        int first_cycle_index = findBehaviorIndex( "exit from cycle phase 0" ); //  4*m; 
        int max_cycle_index = pCD.phenotype.cycle.phases.size();
        if( max_cycle_index > 6 )
        {
            max_cycle_index = 6;
            System.out.println( "Warning: Standardized behaviors only support exit rate from the first 6 phases of a cell cycle!"
                    + "         Ignoring any later phase exit rates." );
        }
        for( int i = 0; i < max_cycle_index; i++ )
        {
            parameters[first_cycle_index + i] = pCD.phenotype.cycle.data.getExitRate( i );
        }

        int apoptosis_index = pCD.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );
        int necrosis_index = pCD.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );

        // apoptosis 
        int apoptosis_param_index = findBehaviorIndex( "apoptosis" );
        parameters[apoptosis_param_index] = pCD.phenotype.death.rates.get( apoptosis_index );

        // necrosis 
        int necrosis_param_index = findBehaviorIndex( "necrosis" );
        parameters[necrosis_param_index] = pCD.phenotype.death.rates.get( necrosis_index );

        // migration speed
        int migration_speed_index = findBehaviorIndex( "migration speed" );
        parameters[migration_speed_index] = pCD.phenotype.motility.migrationSpeed;

        // migration bias 
        int migration_bias_index = findBehaviorIndex( "migration bias" );
        parameters[migration_bias_index] = pCD.phenotype.motility.migrationBias;

        // migration persistence time
        int migration_pt_index = findBehaviorIndex( "migration persistence time" );
        parameters[migration_pt_index] = pCD.phenotype.motility.persistenceTime;

        // chemotactic sensitivities 
        int first_chemotaxis_index = findBehaviorIndex( "chemotactic response to " + microenvironment.densityNames[0] );
        System.arraycopy( pCD.phenotype.motility.chemotacticSensitivities, 0, parameters, first_chemotaxis_index,
                pCD.phenotype.motility.chemotacticSensitivities.length );
        //        std::copy(  pCD.phenotype.motility.chemotactic_sensitivities.begin() ,
        //                    pCD.phenotype.motility.chemotactic_sensitivities.end() ,
        //                    parameters.begin()+first_chemotaxis_index ); 

        // cell-cell adhesion 
        int cca_index = findBehaviorIndex( "cell-cell adhesion" );
        parameters[cca_index] = pCD.phenotype.mechanics.cellCellAdhesionStrength;

        // cell-cell "springs"
        int cca_spring_index = findBehaviorIndex( "cell-cell adhesion elastic constant" );
        parameters[cca_spring_index] = pCD.phenotype.mechanics.attachmentElasticConstant;

        // cell adhesion affinities 
        String search_for1 = "adhesive affinity to " + model.getCellDefinition( 0 ).name;//CellDefinitions_by_type[0].name;
        int first_affinity_index = findBehaviorIndex( search_for1 );
        System.arraycopy( pCD.phenotype.mechanics.cellAdhesionAffinities, 0, parameters, first_affinity_index,
                pCD.phenotype.mechanics.cellAdhesionAffinities.length );
        //        std::copy(  pCD.phenotype.mechanics.cell_adhesion_affinities.begin(), 
        //                    pCD.phenotype.mechanics.cell_adhesion_affinities.end() ,
        //                    parameters.begin()+first_affinity_index ); 

        // max relative maximum adhesion distance 
        int max_adhesion_distance_index = findBehaviorIndex( "relative maximum adhesion distance" );
        parameters[max_adhesion_distance_index] = pCD.phenotype.mechanics.relMaxAdhesionDistance;

        // cell-cell repulsion 
        int ccr_index = findBehaviorIndex( "cell-cell repulsion" );
        parameters[ccr_index] = pCD.phenotype.mechanics.cellCellRepulsionStrength;

        // cell-BM adhesion 
        int cba_index = findBehaviorIndex( "cell-BM adhesion" );
        parameters[cba_index] = pCD.phenotype.mechanics.cellBMAdhesionStrength;

        // cell-BM repulsion 
        int cbr_index = findBehaviorIndex( "cell-BM repulsion" );
        parameters[cbr_index] = pCD.phenotype.mechanics.cellBMRepulsionStrength;

        // dead cell phagocytosis
        int dead_phag_index = findBehaviorIndex( "phagocytose dead cell" );
        parameters[dead_phag_index] = pCD.phenotype.cellInteractions.deadPhagocytosisRate;

        // phagocytosis of each live cell type 
        int first_phagocytosis_index = findBehaviorIndex( "phagocytose " + model.getCellDefinition( 0 ).name );//CellDefinitions_by_type[0].name );
        System.arraycopy( pCD.phenotype.cellInteractions.livePhagocytosisRates, 0, parameters, first_phagocytosis_index,
                pCD.phenotype.cellInteractions.livePhagocytosisRates.length );
        //        std::copy(  pCD.phenotype.cell_interactions.live_phagocytosis_rates.begin(), 
        //                    pCD.phenotype.cell_interactions.live_phagocytosis_rates.end(), 
        //                    parameters.begin()+first_phagocytosis_index );  

        // attack of each live cell type 
        int first_attack_index = findBehaviorIndex( "attack " + model.getCellDefinition( 0 ).name );//CellDefinitions_by_type[0].name );
        System.arraycopy( pCD.phenotype.cellInteractions.attackRates, 0, parameters, first_attack_index,
                pCD.phenotype.cellInteractions.attackRates.length );
        //        std::copy(  pCD.phenotype.cell_interactions.attack_rates.begin(), 
        //                    pCD.phenotype.cell_interactions.attack_rates.end(), 
        //                    parameters.begin()+first_attack_index );    

        // fusion 
        int first_fusion_index = findBehaviorIndex( "fuse to " + model.getCellDefinition( 0 ).name );//CellDefinitions_by_type[0].name );
        System.arraycopy( pCD.phenotype.cellInteractions.fusionRates, 0, parameters, first_fusion_index,
                pCD.phenotype.cellInteractions.fusionRates.length );
        //        std::copy(  pCD.phenotype.cell_interactions.fusion_rates.begin(), 
        //                    pCD.phenotype.cell_interactions.fusion_rates.end(), 
        //                    parameters.begin()+first_fusion_index );    

        // transformation 
        int first_transformation_index = findBehaviorIndex( "transform to " + model.getCellDefinition( 0 ).name );
        System.arraycopy( pCD.phenotype.cellTransformations.transformationRates, 0, parameters, first_transformation_index,
                pCD.phenotype.cellTransformations.transformationRates.length );
        //        std::copy(  pCD.phenotype.cell_transformations.transformation_rates.begin(), 
        //                    pCD.phenotype.cell_transformations.transformation_rates.end(), 
        //                    parameters.begin()+first_transformation_index );    

        // custom behavior
        int first_custom_ind = findBehaviorIndex( "custom 0" );
        if( first_custom_ind >= 0 )
        {
            for( int nc = 0; nc < pCell.customData.variables.size(); nc++ )
            {
                parameters[first_custom_ind + nc] = pCD.custom_data.variables.get( nc ).value;
            }
        }

        // is the cell movable / not movable 
        int movable_ind = findBehaviorIndex( "is_movable" );
        if( pCD.isMovable == true )
        {
            parameters[movable_ind] = 1;
        }
        else
        {
            parameters[movable_ind] = 0;
        }

        // vector of immunogenicity behaviors 
        int start_immunogenicity_ind = findBehaviorIndex( "immunogenicity to " + model.getCellDefinition( 0 ).name );//[0].name );
        System.arraycopy( pCD.phenotype.cellInteractions.immunogenicities, 0, parameters, start_immunogenicity_ind,
                pCD.phenotype.cellInteractions.immunogenicities.length );
        //        std::copy( pCD.phenotype.cell_interactions.immunogenicities.begin(),
        //                   pCD.phenotype.cell_interactions.immunogenicities.end(), 
        //                   parameters.begin()+start_immunogenicity_ind );  


        // set cell attachment rate  
        int attachment_rate_ind = findBehaviorIndex( "cell attachment rate" );
        parameters[attachment_rate_ind] = pCD.phenotype.mechanics.attachmentRate;

        // set cell detachment rate  
        int detachment_rate_ind = findBehaviorIndex( "cell detachment rate" );
        parameters[detachment_rate_ind] = pCD.phenotype.mechanics.detachmentRate;

        // maximum number of cell attachments 
        int max_attachments_ind = findBehaviorIndex( "maximum number of cell attachments" );
        parameters[max_attachments_ind] = pCD.phenotype.mechanics.maxAttachments;

        // cell damage rate (effector attack)
        int damage_rate_ind = findBehaviorIndex( "damage rate" );
        parameters[damage_rate_ind] = pCD.phenotype.cellInteractions.damageRate;

        return parameters;
    }

    public double getSingleBaseBehavior(Cell pCell, int index)
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();
        Model model = pCell.getModel();
        int m = microenvironment.numberDensities();
        int n = model.getDefinitionsCount();

        CellDefinition pCD = model.getCellDefinition( pCell.typeName );

        if( index < 0 )
        {
            System.out.println(
                    "Warning: attempted to get behavior with unknown index " + "         I'm ignoring it, but you should fix it!" );
            return 0.0;
        }

        // substrate-related behaviors 

        // first m entries are secretion 
        int first_secretion_index = findBehaviorIndex( microenvironment.densityNames[0] + " secretion" ); // 0; 
        if( index >= first_secretion_index && index < first_secretion_index + m )
        {
            return pCD.phenotype.secretion.secretionRates[index - first_secretion_index];
        }

        // next m entries are secretion targets
        int first_secretion_target_index = findBehaviorIndex( microenvironment.densityNames[0] + " secretion target" ); // m; 
        if( index >= first_secretion_target_index && index < first_secretion_target_index + m )
        {
            return pCD.phenotype.secretion.saturationDensities[index - first_secretion_target_index];
        }

        // next m entries are uptake rates
        int first_uptake_index = findBehaviorIndex( microenvironment.densityNames[0] + " uptake" ); // 2*m; 
        if( index >= first_uptake_index && index < first_uptake_index + m )
        {
            return pCD.phenotype.secretion.uptakeRates[index - first_uptake_index];
        }

        // next m entries are net export rates 
        int first_export_index = findBehaviorIndex( microenvironment.densityNames[0] + " export" ); //  3*m; 
        if( index >= first_export_index && index < first_export_index + m )
        {
            return pCD.phenotype.secretion.netExportRates[index - first_export_index];
        }

        // cycle entry (exit from phase 0) and exit from up to 5 more phases 
        int first_cycle_index = findBehaviorIndex( "exit from cycle phase 0" ); //  4*m; 
        int max_cycle_index = pCD.phenotype.cycle.phases.size();
        if( max_cycle_index > 6 )
        {
            max_cycle_index = 6;
            System.out.println( "Warning: Standardized behaviors only support exit rate from the first 6 phases of a cell cycle!"
                    + "         Ignoring any later phase exit rates." );
        }
        if( index >= first_cycle_index && index < first_cycle_index + 6 )
        {
            int ind = index - first_cycle_index;
            if( ind < max_cycle_index )
            {
                return pCD.phenotype.cycle.data.getExitRate( ind );
            }
            return 0.0;
        }

        int apoptosis_index = pCD.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );
        int necrosis_index = pCD.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );

        int apop_param_index = findBehaviorIndex( "apoptosis" );
        int necr_param_index = findBehaviorIndex( "necrosis" );

        // apoptosis 
        if( index == apop_param_index )
        {
            return pCD.phenotype.death.rates.get( apoptosis_index );
        }

        // necrosis 
        if( index == necr_param_index )
        {
            return pCD.phenotype.death.rates.get( necrosis_index );
        }

        // migration speed
        int migr_spd_index = findBehaviorIndex( "migration speed" );
        if( index == migr_spd_index )
        {
            return pCD.phenotype.motility.migrationSpeed;
        }

        // migration bias 
        int migr_bias_index = findBehaviorIndex( "migration bias" );
        if( index == migr_bias_index )
        {
            return pCD.phenotype.motility.migrationBias;
        }

        // migration persistence time
        int migr_pt_index = findBehaviorIndex( "migration persistence time" );
        if( index == migr_pt_index )
        {
            return pCD.phenotype.motility.persistenceTime;
        }

        // chemotactic sensitivities 
        int first_chemotaxis_index = findBehaviorIndex( "chemotactic response to " + microenvironment.densityNames[0] );
        if( index >= first_chemotaxis_index && index < first_chemotaxis_index + m )
        {
            return pCD.phenotype.motility.chemotacticSensitivities[index - first_chemotaxis_index];
        }

        // cell-cell adhesion 
        int cca_index = findBehaviorIndex( "cell-cell adhesion" );
        if( index == cca_index )
        {
            return pCD.phenotype.mechanics.cellCellAdhesionStrength;
        }

        // cell-cell "springs"
        int cca_spring_index = findBehaviorIndex( "cell-cell adhesion elastic constant" );
        if( index == cca_spring_index )
        {
            return pCD.phenotype.mechanics.attachmentElasticConstant;
        }

        // cell adhesion affinities 
        int first_affinity_index = findBehaviorIndex( "adhesive affinity to " + model.getCellDefinition( 0 ).name );//  CellDefinitions_by_type[0].name );
        if( index >= first_affinity_index && index < first_affinity_index + n )
        {
            return pCD.phenotype.mechanics.cellAdhesionAffinities[index - first_affinity_index];
        }

        // max relative maximum adhesion distance 
        int max_adh_index = findBehaviorIndex( "relative maximum adhesion distance" );
        if( index == max_adh_index )
        {
            return pCD.phenotype.mechanics.relMaxAdhesionDistance;
        }

        // cell-cell repulsion 
        int ccr_index = findBehaviorIndex( "cell-cell repulsion" );
        if( index == ccr_index )
        {
            return pCD.phenotype.mechanics.cellCellRepulsionStrength;
        }

        // cell-BM adhesion 
        int cba_index = findBehaviorIndex( "cell-BM adhesion" );
        if( index == cba_index )
        {
            return pCD.phenotype.mechanics.cellBMAdhesionStrength;
        }

        // cell-BM repulsion 
        int cbr_index = findBehaviorIndex( "cell-BM repulsion" );
        if( index == cbr_index )
        {
            return pCD.phenotype.mechanics.cellBMRepulsionStrength;
        }

        // dead cell phagocytosis
        int dead_phag_index = findBehaviorIndex( "phagocytose dead cell" );
        if( index == dead_phag_index )
        {
            return pCD.phenotype.cellInteractions.deadPhagocytosisRate;
        }

        // phagocytosis of each live cell type 
        int first_phagocytosis_index = findBehaviorIndex( "phagocytose " + model.getCellDefinition( 0 ).name );//CellDefinitions_by_type[0].name );
        if( index >= first_phagocytosis_index && index < first_phagocytosis_index + n )
        {
            return pCD.phenotype.cellInteractions.livePhagocytosisRates[index - first_phagocytosis_index];
        }

        // attack of each live cell type 
        int first_attack_index = findBehaviorIndex( "attack " + model.getCellDefinition( 0 ).name );//CellDefinitions_by_type[0].name );
        if( index >= first_attack_index && index < first_attack_index + n )
        {
            return pCD.phenotype.cellInteractions.attackRates[index - first_attack_index];
        }

        // fusion 
        int first_fusion_index = findBehaviorIndex( "fuse to " + model.getCellDefinition( 0 ).name );//CellDefinitions_by_type[0].name );
        if( index >= first_fusion_index && index < first_fusion_index + n )
        {
            return pCD.phenotype.cellInteractions.fusionRates[index - first_fusion_index];
        }

        // transformation 
        int first_transformation_index = findBehaviorIndex( "transform to " + model.getCellDefinition( 0 ).name );//CellDefinitions_by_type[0].name );
        if( index >= first_transformation_index && index < first_transformation_index + n )
        {
            return pCD.phenotype.cellTransformations.transformationRates[index - first_transformation_index];
        }

        // custom behavior
        int first_custom_ind = findBehaviorIndex( "custom 0" );
        int max_custom_ind = first_custom_ind + pCell.customData.variables.size();
        if( first_custom_ind >= 0 && index >= first_custom_ind && index < max_custom_ind )
        {
            return pCD.custom_data.variables.get( index - first_custom_ind ).value;
        }

        // is the cell movable / not movable 
        int movable_ind = findBehaviorIndex( "is_movable" );
        if( index == movable_ind )
        {
            if( pCD.isMovable == true )
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        }

        // vector of immunogenicity behaviors 
        int start_immunogenicity_ind = findBehaviorIndex( "immunogenicity to " + model.getCellDefinition( 0 ).name );//CellDefinitions_by_type[0].name );
        int max_immunogenicity_ind = start_immunogenicity_ind + n;
        if( start_immunogenicity_ind > -1 && index >= start_immunogenicity_ind && index < max_immunogenicity_ind )
        {
            return pCD.phenotype.cellInteractions.immunogenicities[index - start_immunogenicity_ind];
        }


        // set cell attachment rate  
        int attachment_rate_ind = findBehaviorIndex( "cell attachment rate" );
        if( index == attachment_rate_ind )
        {
            return pCD.phenotype.mechanics.attachmentRate;
        }

        // set cell detachment rate  
        int detachment_rate_ind = findBehaviorIndex( "cell detachment rate" );
        if( index == detachment_rate_ind )
        {
            return pCD.phenotype.mechanics.detachmentRate;
        }

        // maximum number of cell attachments 
        int max_attachments_ind = findBehaviorIndex( "maximum number of cell attachments" );
        if( index == max_attachments_ind )
        {
            return pCD.phenotype.mechanics.maxAttachments;
        }

        // cell damage rate (effector attack)
        int damage_rate_ind = findBehaviorIndex( "damage rate" );
        if( index == damage_rate_ind )
        {
            return pCD.phenotype.cellInteractions.damageRate;
        }

        return -1;
    }

    public double getSingleBaseBehavior(Model model, CellDefinition pCD, String name)
    {
        return getSingleBaseBehavior( model, pCD, findBehaviorIndex( name ) );
    }

    public double getSingleBaseBehavior(Model model, CellDefinition pCD, int index)
    {
        Microenvironment microenvironment = pCD.getMicroenvironment();
        int m = microenvironment.numberDensities();
        int n = model.getDefinitionsCount();

        // CellDefinition* pCD = find_CellDefinition( pCell.type_name );     

        if( index < 0 )
        {
            System.out.println(
                    "Warning: attempted to get behavior with unknown index " + "         I'm ignoring it, but you should fix it!" );
            return 0.0;
        }

        // substrate-related behaviors 

        // first m entries are secretion 
        int first_secretion_index = findBehaviorIndex( microenvironment.densityNames[0] + " secretion" ); // 0; 
        if( index >= first_secretion_index && index < first_secretion_index + m )
        {
            return pCD.phenotype.secretion.secretionRates[index - first_secretion_index];
        }

        // next m entries are secretion targets
        int first_secretion_target_index = findBehaviorIndex( microenvironment.densityNames[0] + " secretion target" ); // m; 
        if( index >= first_secretion_target_index && index < first_secretion_target_index + m )
        {
            return pCD.phenotype.secretion.saturationDensities[index - first_secretion_target_index];
        }

        // next m entries are uptake rates
        int first_uptake_index = findBehaviorIndex( microenvironment.densityNames[0] + " uptake" ); // 2*m; 
        if( index >= first_uptake_index && index < first_uptake_index + m )
        {
            return pCD.phenotype.secretion.uptakeRates[index - first_uptake_index];
        }

        // next m entries are net export rates 
        int first_export_index = findBehaviorIndex( microenvironment.densityNames[0] + " export" ); //  3*m; 
        if( index >= first_export_index && index < first_export_index + m )
        {
            return pCD.phenotype.secretion.netExportRates[index - first_export_index];
        }

        // cycle entry (exit from phase 0) and exit from up to 5 more phases 
        int first_cycle_index = findBehaviorIndex( "exit from cycle phase 0" ); //  4*m; 
        int max_cycle_index = pCD.phenotype.cycle.phases.size();
        if( max_cycle_index > 6 )
        {
            max_cycle_index = 6;
            System.out.println( "Warning: Standardized behaviors only support exit rate from the first 6 phases of a cell cycle!\n"
                    + "         Ignoring any later phase exit rates." );
        }
        if( index >= first_cycle_index && index < first_cycle_index + 6 )
        {
            int ind = index - first_cycle_index;
            if( ind < max_cycle_index )
            {
                return pCD.phenotype.cycle.data.getExitRate( ind );
            }
            return 0.0;
        }

        int apoptosis_index = pCD.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );
        int necrosis_index = pCD.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );

        int apop_param_index = findBehaviorIndex( "apoptosis" );
        int necr_param_index = findBehaviorIndex( "necrosis" );

        // apoptosis 
        if( index == apop_param_index )
        {
            return pCD.phenotype.death.rates.get( apoptosis_index );
        }

        // necrosis 
        if( index == necr_param_index )
        {
            return pCD.phenotype.death.rates.get( necrosis_index );
        }

        // migration speed
        int migr_spd_index = findBehaviorIndex( "migration speed" );
        if( index == migr_spd_index )
        {
            return pCD.phenotype.motility.migrationSpeed;
        }

        // migration bias 
        int migr_bias_index = findBehaviorIndex( "migration bias" );
        if( index == migr_bias_index )
        {
            return pCD.phenotype.motility.migrationBias;
        }

        // migration persistence time
        int migr_pt_index = findBehaviorIndex( "migration persistence time" );
        if( index == migr_pt_index )
        {
            return pCD.phenotype.motility.persistenceTime;
        }

        // chemotactic sensitivities 
        int first_chemotaxis_index = findBehaviorIndex( "chemotactic response to " + microenvironment.densityNames[0] );
        if( index >= first_chemotaxis_index && index < first_chemotaxis_index + m )
        {
            return pCD.phenotype.motility.chemotacticSensitivities[index - first_chemotaxis_index];
        }

        // cell-cell adhesion 
        int cca_index = findBehaviorIndex( "cell-cell adhesion" );
        if( index == cca_index )
        {
            return pCD.phenotype.mechanics.cellCellAdhesionStrength;
        }

        // cell-cell "springs"
        int cca_spring_index = findBehaviorIndex( "cell-cell adhesion elastic constant" );
        if( index == cca_spring_index )
        {
            return pCD.phenotype.mechanics.attachmentElasticConstant;
        }

        // cell adhesion affinities 
        int first_affinity_index = findBehaviorIndex( "adhesive affinity to " + model.getCellDefinition( 0 ).name );
        if( index >= first_affinity_index && index < first_affinity_index + n )
        {
            return pCD.phenotype.mechanics.cellAdhesionAffinities[index - first_affinity_index];
        }

        // max relative maximum adhesion distance 
        int max_adh_index = findBehaviorIndex( "relative maximum adhesion distance" );
        if( index == max_adh_index )
        {
            return pCD.phenotype.mechanics.relMaxAdhesionDistance;
        }

        // cell-cell repulsion 
        int ccr_index = findBehaviorIndex( "cell-cell repulsion" );
        if( index == ccr_index )
        {
            return pCD.phenotype.mechanics.cellCellRepulsionStrength;
        }

        // cell-BM adhesion 
        int cba_index = findBehaviorIndex( "cell-BM adhesion" );
        if( index == cba_index )
        {
            return pCD.phenotype.mechanics.cellBMAdhesionStrength;
        }

        // cell-BM repulsion 
        int cbr_index = findBehaviorIndex( "cell-BM repulsion" );
        if( index == cbr_index )
        {
            return pCD.phenotype.mechanics.cellBMRepulsionStrength;
        }

        // dead cell phagocytosis
        int dead_phag_index = findBehaviorIndex( "phagocytose dead cell" );
        if( index == dead_phag_index )
        {
            return pCD.phenotype.cellInteractions.deadPhagocytosisRate;
        }

        // phagocytosis of each live cell type 
        int first_phagocytosis_index = findBehaviorIndex( "phagocytose " + model.getCellDefinition( 0 ).name );
        if( index >= first_phagocytosis_index && index < first_phagocytosis_index + n )
        {
            return pCD.phenotype.cellInteractions.livePhagocytosisRates[index - first_phagocytosis_index];
        }

        // attack of each live cell type 
        int first_attack_index = findBehaviorIndex( "attack " + model.getCellDefinition( 0 ).name );
        if( index >= first_attack_index && index < first_attack_index + n )
        {
            return pCD.phenotype.cellInteractions.attackRates[index - first_attack_index];
        }

        // fusion 
        int first_fusion_index = findBehaviorIndex( "fuse to " + model.getCellDefinition( 0 ).name );
        if( index >= first_fusion_index && index < first_fusion_index + n )
        {
            return pCD.phenotype.cellInteractions.fusionRates[index - first_fusion_index];
        }

        // transformation 
        int first_transformation_index = findBehaviorIndex( "transform to " + model.getCellDefinition( 0 ).name );
        if( index >= first_transformation_index && index < first_transformation_index + n )
        {
            return pCD.phenotype.cellTransformations.transformationRates[index - first_transformation_index];
        }

        // custom behavior
        int first_custom_ind = findBehaviorIndex( "custom 0" );
        int max_custom_ind = first_custom_ind + pCD.custom_data.variables.size();
        if( first_custom_ind >= 0 && index >= first_custom_ind && index < max_custom_ind )
        {
            return pCD.custom_data.variables.get( index - first_custom_ind ).value;
        }

        // is the cell movable / not movable 
        int movable_ind = findBehaviorIndex( "is_movable" );
        if( index == movable_ind )
        {
            if( pCD.isMovable == true )
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        }

        // vector of immunogenicity behaviors 
        int start_immunogenicity_ind = findBehaviorIndex( "immunogenicity to " + model.getCellDefinition( 0 ).name );
        int max_immunogenicity_ind = start_immunogenicity_ind + n;
        if( start_immunogenicity_ind > -1 && index >= start_immunogenicity_ind && index < max_immunogenicity_ind )
        {
            return pCD.phenotype.cellInteractions.immunogenicities[index - start_immunogenicity_ind];
        }


        // set cell attachment rate  
        int attachment_rate_ind = findBehaviorIndex( "cell attachment rate" );
        if( index == attachment_rate_ind )
        {
            return pCD.phenotype.mechanics.attachmentRate;
        }

        // set cell detachment rate  
        int detachment_rate_ind = findBehaviorIndex( "cell detachment rate" );
        if( index == detachment_rate_ind )
        {
            return pCD.phenotype.mechanics.detachmentRate;
        }

        // maximum number of cell attachments 
        int max_attachments_ind = findBehaviorIndex( "maximum number of cell attachments" );
        if( index == max_attachments_ind )
        {
            return pCD.phenotype.mechanics.maxAttachments;
        }

        // cell damage rate (effector attack)
        int damage_rate_ind = findBehaviorIndex( "damage rate" );
        if( index == damage_rate_ind )
        {
            return pCD.phenotype.cellInteractions.damageRate;
        }

        return -1;
    }

    public double get_single_base_behavior(Cell pCell, String name)
    {
        return getSingleBaseBehavior( pCell, findBehaviorIndex( name ) );
    }

    public double[] get_base_behaviors(Cell pCell, int[] indices)
    {
        double[] parameters = new double[indices.length];//( indices.size() , 0.0 ); 
        for( int n = 0; n < indices.length; n++ )
        {
            parameters[n] = getSingleBaseBehavior( pCell, indices[n] );
        }
        return parameters;
    }

    public double[] get_base_behaviors(Cell pCell, String[] names)
    {
        double[] parameters = new double[names.length];//;( names.length , 0.0 ); 
        for( int n = 0; n < names.length; n++ )
        {
            parameters[n] = get_single_base_behavior( pCell, names[n] );
        }
        return parameters;
    }
}
