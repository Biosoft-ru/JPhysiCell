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

    //signal indexes
    private int startSubstrate;
    private int startSubstrateIntra;
    private int startSubstrateGrad;
    private int pressure;
    private int volume;
    private int contact;
    private int contactLive;
    private int contactDead;
    private int contactBM;
    private int damage;
    private int dead;
    private int totalAttackTime;
    private int time;
    private int startCustom;
    private int apoptotic;
    private int necrotic;

    //behavior indexes
    private int first_secretion_index;
    private int first_secretion_target_index;
    private int first_uptake_index;
    private int first_export_index;
    private int first_cycle_index;
    private int apoptosis_parameter_index;
    private int necrosis_parameter_index;
    private int migration_speed_index;
    private int migration_bias_index;
    private int persistence_time_index;
    private int first_chemotaxis_index;
    private int cca_index;
    private int elastic_index;
    private int first_affinity_index;
    private int max_adh_distance_index;
    private int ccr_index;
    private int cba_index;
    private int cbr_index;
    private int dead_phago_index;
    private int first_phagocytosis_index;
    private int first_attack_index;
    private int first_fusion_index;
    private int first_transformation_index;
    private int first_custom_ind;
    private int movable_ind;
    private int first_immunogenicity_index;
    private int attachment_rate_ind;
    private int detachment_rate_ind;
    private int max_attachments_ind;
    private int damage_rate_ind;

    private int signalsNumber;
    private int behaviorsNumber;

    public double getSingleSignal(Cell cell, int index) throws IllegalArgumentException
    {
        Microenvironment microenvironment = cell.getMicroenvironment();
        int m = microenvironment.numberDensities();
        int n = cell.getModel().getDefinitionsCount();

        if( startSubstrate <= index && index < startSubstrate + m )
            return scale( cell.nearest_density_vector()[index - startSubstrate], index );
        else if( startSubstrateIntra <= index && index < startSubstrateIntra + m )
            return scale( cell.phenotype.molecular.internSubstrates[index - startSubstrateIntra] / cell.phenotype.volume.total, index );
        else if( startSubstrateGrad <= index && index < startSubstrateGrad + m )
            return scale( VectorUtil.norm( cell.nearest_gradient( index - startSubstrateGrad ) ), index );
        else if( index == pressure )
            return scale( cell.state.simplePressure, index );
        else if( index == volume )
            return scale( cell.phenotype.volume.total, index );
        else if( contact <= index && index < contact + n + 2 )
        {
            int[] counts = new int[n];
            int dead_cells = 0;
            int live_cells = 0;
            for( Cell pC : cell.state.neighbors ) // process all neighbors 
            {
                if( pC.phenotype.death.dead )
                    dead_cells++;
                else
                    live_cells++;
                counts[pC.type] += 1;
            }

            if( index < contact + n )
                return scale( counts[index - contact], index );
            else if( index == contactLive )
                return scale( live_cells, index );
            else if( index == contactDead )
                return scale( dead_cells, index );
        }
        else if( index == contactBM )
            return scale( digitize( cell.state.contactWithBasementMembrane ), index );
        else if( index == damage )
            return scale( cell.state.damage, index );
        else if( index == dead )
            return scale( digitize( cell.phenotype.death.dead ), index );
        else if( index == totalAttackTime )
            return scale( cell.state.totalAttackTime, index );
        else if( index == time )
            return scale( cell.getMicroenvironment().time, index );
        else if( startCustom > -1 && index >= startCustom && index < startCustom + cell.customData.variables.size() )
            return scale( cell.customData.variables.get( index - startCustom ).value, index );
        else if( index == apoptotic )
            return digitize( cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.apoptotic );
        else if( index == necrotic )
            return digitize( cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                    || cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                    || cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic );

        throw new IllegalArgumentException( "Incorrect signal with index " + index );
    }

    private double scale(double val, int index)
    {
        return val / signalScales[index];
    }

    private double digitize(boolean val)
    {
        return val ? 1 : 0;
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

    public int findSignalIndex(String name)
    {
        if( signalToIndex.containsKey( name ) )
            return signalToIndex.get( name );
        throw new IllegalArgumentException( "having trouble finding " + name );
    }

    public double getSingleSignal(Cell pCell, String name)
    {
        return getSingleSignal( pCell, findSignalIndex( name ) );
    }

    private void registerSignal(String ... names)
    {
        indexToSignal.put( signalsNumber, names[0] );
        for( String name : names )
            signalToIndex.put( name, signalsNumber );

        signalsNumber++;
    }

    private void clearSignals()
    {
        signalToIndex.clear();
        indexToSignal.clear();
        signalsNumber = 0;
    }

    private void clearBehaviors()
    {
        behaviorToInt.clear();
        indexToBehavior.clear();
        behaviorsNumber = 0;
    }

    private void registerBehavior(String ... names)
    {
        indexToBehavior.put( behaviorsNumber, names[0] );
        for( String name : names )
            behaviorToInt.put( name, behaviorsNumber );
        behaviorsNumber++;
    }

    public void setupDictionaries(Model model)
    {
        clearSignals();

        Microenvironment m = model.getMicroenvironment();
        int densNumber = m.numberDensities();

        for( int i = 0; i < densNumber; i++ )
            registerSignal( m.densityNames[i] );

        for( int i = 0; i < densNumber; i++ )
            registerSignal( "intracellular " + m.densityNames[i], "internalized " + m.densityNames[i] );

        for( int i = 0; i < densNumber; i++ )
            registerSignal( m.densityNames[i] + " gradient", "grad(" + m.densityNames[i] + ")", "gradient of " + m.densityNames[i] );

        registerSignal( "pressure" );
        registerSignal( "volume" );

        for( CellDefinition cd : model.getCellDefinitions() )
            registerSignal( "contact with cell type " + cd.type, "contact with " + cd.name );

        registerSignal( "contact with live cell", "contact with live cells" );
        registerSignal( "contact with dead cell", "contact with dead cells" );
        registerSignal( "contact with basement membrane", "contact with BM" );
        registerSignal( "damage" );
        registerSignal( "dead", "is dead" );
        registerSignal( "total attack time" );
        registerSignal( "time", "current time", "global time" );

        CellDefinition def = model.getCellDefinition( 0 );
        for( int i = 0; i < def.custom_data.variables.size(); i++ )
        {
            String varName = def.custom_data.variables.get( i ).name;
            registerSignal( "custom:" + varName, "custom: " + varName, "custom " + i );
        }

        registerSignal( "apoptotic", "is_apoptotic" );
        registerSignal( "necrotic", "is_necrotic" );


        startSubstrate = findSignalIndex( m.densityNames[0] );
        startSubstrateIntra = findSignalIndex( "intracellular " + m.densityNames[0] );
        startSubstrateGrad = findSignalIndex( m.densityNames[0] + " gradient" );
        pressure = findSignalIndex( "pressure" );
        volume = findSignalIndex( "volume" );
        contact = findSignalIndex( "contact with " + def.name );
        contactLive = findSignalIndex( "contact with live cell" );
        contactDead = findSignalIndex( "contact with dead cell" );
        contactBM = findSignalIndex( "contact with basement membrane" );
        damage = findSignalIndex( "damage" );
        dead = findSignalIndex( "dead" );
        totalAttackTime = findSignalIndex( "total attack time" );
        time = findSignalIndex( "time" );
        startCustom = findSignalIndex( "custom:" + def.custom_data.variables.get( 0 ).name );
        apoptotic = findSignalIndex( "apoptotic" );
        necrotic = findSignalIndex( "necrotic" );

        clearBehaviors();

        for( int i = 0; i < densNumber; i++ )
            registerBehavior( m.densityNames[i] + " secretion" );

        for( int i = 0; i < densNumber; i++ )
            registerBehavior( m.densityNames[i] + " secretion target", m.densityNames[i] + " secretion saturation density" );

        for( int i = 0; i < densNumber; i++ )
            registerBehavior( m.densityNames[i] + " uptake" );

        for( int i = 0; i < densNumber; i++ )
            registerBehavior( m.densityNames[i] + " export" );

        registerBehavior( "cycle entry", "exit from cycle phase 0" );

        for( int i = 1; i < 6; i++ )
            registerBehavior( "exit from cycle phase " + i );

        registerBehavior( "apoptosis" );
        registerBehavior( "necrosis" );
        registerBehavior( "migration speed" );
        registerBehavior( "migration bias" );
        registerBehavior( "migration persistence time" );

        for( int i = 0; i < densNumber; i++ )
            registerBehavior( "chemotactic response to " + m.densityNames[i], "chemotactic sensitivity to " + m.densityNames[i] );

        registerBehavior( "cell-cell adhesion" );
        registerBehavior( "cell-cell adhesion elastic constant" );

        for( CellDefinition cd : model.getCellDefinitions() )
            registerBehavior( "adhesive affinity to " + cd.name, "adhesive affinity to cell type " + cd.type );

        registerBehavior( "relative maximum adhesion distance" );
        registerBehavior( "cell-cell repulsion" );
        registerBehavior( "cell-BM adhesion", "cell-membrane adhesion" );
        registerBehavior( "cell-BM repulsion", "cell-membrane repulsion" );
        registerBehavior( "phagocytose dead cell", "phagocytosis of dead cell", "phagocytosis of dead cells" );

        for( CellDefinition cd : model.getCellDefinitions() )
            registerBehavior( "phagocytose " + cd.name, "phagocytose cell type " + cd.type, "phagocytosis of " + cd.type );

        for( CellDefinition cd : model.getCellDefinitions() )
            registerBehavior( "attack " + cd.name, "attack cell type " + cd.type );

        for( CellDefinition cd : model.getCellDefinitions() )
            registerBehavior( "fuse to " + cd.name, "fuse to cell type " + cd.type );

        for( CellDefinition cd : model.getCellDefinitions() )
            registerBehavior( "transform to " + cd.name, "transform to cell type " + cd.type );

        for( int i = 0; i < def.custom_data.variables.size(); i++ )
        {
            String varName = def.custom_data.variables.get( i ).getName();
            registerBehavior( "custom:" + varName, "custom: " + varName, "custom " + i );
        }

        registerBehavior( "is_movable", "movable", "is movable" );

        for( CellDefinition cd : model.getCellDefinitions() )
            registerBehavior( "immunogenicity to " + cd.name, "immunogenicity to cell type " + cd.type );

        registerBehavior( "cell attachment rate" );
        registerBehavior( "cell detachment rate" );
        registerBehavior( "maximum number of cell attachments" );
        registerBehavior( "damage rate" );

        signalScales = VectorUtil.resize( signalScales, indexToSignal.size(), 1.0 );

        first_secretion_index = findBehaviorIndex( m.densityNames[0] + " secretion" );
        first_secretion_target_index = findBehaviorIndex( m.densityNames[0] + " secretion target" );
        first_uptake_index = findBehaviorIndex( m.densityNames[0] + " uptake" );
        first_export_index = findBehaviorIndex( m.densityNames[0] + " export" );
        first_cycle_index = findBehaviorIndex( "cycle entry" );
        apoptosis_parameter_index = findBehaviorIndex( "apoptosis" );
        necrosis_parameter_index = findBehaviorIndex( "necrosis" );
        migration_speed_index = findBehaviorIndex( "migration speed" );
        migration_bias_index = findBehaviorIndex( "migration bias" );
        persistence_time_index = findBehaviorIndex( "migration persistence time" );
        first_chemotaxis_index = findBehaviorIndex( "chemotactic response to " + m.densityNames[0] );
        cca_index = findBehaviorIndex( "cell-cell adhesion" );
        elastic_index = findBehaviorIndex( "cell-cell adhesion elastic constant" );
        first_affinity_index = findBehaviorIndex( "adhesive affinity to " + def.name );
        max_adh_distance_index = findBehaviorIndex( "relative maximum adhesion distance" );
        ccr_index = findBehaviorIndex( "cell-cell repulsion" );
        cba_index = findBehaviorIndex( "cell-BM adhesion" );
        cbr_index = findBehaviorIndex( "cell-BM repulsion" );
        dead_phago_index = findBehaviorIndex( "phagocytose dead cell" );
        first_phagocytosis_index = findBehaviorIndex( "phagocytose " + def.name );
        first_attack_index = findBehaviorIndex( "attack " + def.name );
        first_fusion_index = findBehaviorIndex( "fuse to " + def.name );
        first_transformation_index = findBehaviorIndex( "transform to " + def.name );
        first_custom_ind = findBehaviorIndex( "custom:" + def.custom_data.variables.get( 0 ).name );
        movable_ind = findBehaviorIndex( "is_movable" );
        first_immunogenicity_index = findBehaviorIndex( "immunogenicity to " + def.name );
        attachment_rate_ind = findBehaviorIndex( "cell attachment rate" );
        detachment_rate_ind = findBehaviorIndex( "cell detachment rate" );
        max_attachments_ind = findBehaviorIndex( "maximum number of cell attachments" );
        damage_rate_ind = findBehaviorIndex( "damage rate" );
    }

    void setSingleBehavior(Cell cell, int index, double parameter) throws Exception
    {
        Microenvironment microenvironment = cell.getMicroenvironment();
        Model model = cell.getModel();
        int m = microenvironment.numberDensities();
        int n = model.getDefinitionsCount();

        if( index >= first_secretion_index && index < first_secretion_index + m )
            cell.phenotype.secretion.secretionRates[index - first_secretion_index] = parameter;
        else if( index >= first_secretion_target_index && index < first_secretion_target_index + m )
            cell.phenotype.secretion.saturationDensities[index - first_secretion_target_index] = parameter;
        else if( index >= first_uptake_index && index < first_uptake_index + m )
            cell.phenotype.secretion.uptakeRates[index - first_uptake_index] = parameter;
        else if( index >= first_export_index && index < first_export_index + m )
            cell.phenotype.secretion.netExportRates[index - first_export_index] = parameter;
        else if( index >= first_cycle_index && index < first_cycle_index + 6 )
        {
            int max_cycle_index = cell.phenotype.cycle.phases.size();
            if( index < first_cycle_index + max_cycle_index )
                cell.phenotype.cycle.data.setExitRate( index - first_cycle_index, parameter );
            else
                throw new Exception( "Warning: Attempted to set a cycle exit rate outside the bounds of the cell's cycle model" );
        }

        else if( index == apoptosis_parameter_index )
            cell.phenotype.death.rates.set( cell.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model ),
                    parameter );
        else if( index == necrosis_parameter_index )
            cell.phenotype.death.rates.set( cell.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model ),
                    parameter );
        else if( index == migration_speed_index )
            cell.phenotype.motility.migrationSpeed = parameter;
        else if( index == migration_bias_index )
            cell.phenotype.motility.migrationBias = parameter;
        else if( index == persistence_time_index )
            cell.phenotype.motility.persistenceTime = parameter;
        else if( index >= first_chemotaxis_index && index < first_chemotaxis_index + m )
            cell.phenotype.motility.chemotacticSensitivities[index - first_chemotaxis_index] = parameter;
        else if( index == cca_index )
            cell.phenotype.mechanics.cellCellAdhesionStrength = parameter;
        else if( index == elastic_index )
            cell.phenotype.mechanics.attachmentElasticConstant = parameter;
        else if( index >= first_affinity_index && index < first_affinity_index + n )
            cell.phenotype.mechanics.cellAdhesionAffinities[index - first_affinity_index] = parameter;
        else if( index == max_adh_distance_index )
            cell.phenotype.mechanics.relMaxAdhesionDistance = parameter;
        else if( index == ccr_index )
            cell.phenotype.mechanics.cellCellRepulsionStrength = parameter;
        else if( index == cba_index )
            cell.phenotype.mechanics.cellBMAdhesionStrength = parameter;
        else if( index == cbr_index )
            cell.phenotype.mechanics.cellBMRepulsionStrength = parameter;
        else if( index == dead_phago_index )
            cell.phenotype.cellInteractions.deadPhagocytosisRate = parameter;
        else if( index >= first_phagocytosis_index && index < first_phagocytosis_index + n )
            cell.phenotype.cellInteractions.livePhagocytosisRates[index - first_phagocytosis_index] = parameter;
        else if( index >= first_attack_index && index < first_attack_index + n )
            cell.phenotype.cellInteractions.attackRates[index - first_attack_index] = parameter;
        else if( index >= first_fusion_index && index < first_fusion_index + n )
            cell.phenotype.cellInteractions.fusionRates[index - first_fusion_index] = parameter;
        else if( index >= first_transformation_index && index < first_transformation_index + n )
            cell.phenotype.cellTransformations.transformationRates[index - first_transformation_index] = parameter;
        else if( first_custom_ind >= 0 && index >= first_custom_ind && index < first_custom_ind + cell.customData.variables.size() )
            cell.customData.variables.get( index - first_custom_ind ).value = parameter;
        else if( index == movable_ind )
            cell.isMovable = parameter > 0.5;
        else if( index >= first_immunogenicity_index && index < first_immunogenicity_index + n )
            cell.phenotype.cellInteractions.immunogenicities[index - first_immunogenicity_index] = parameter;
        else if( index == attachment_rate_ind )
            cell.phenotype.mechanics.attachmentRate = parameter;
        else if( index == detachment_rate_ind )
            cell.phenotype.mechanics.detachmentRate = parameter;
        else if( index == max_attachments_ind )
            cell.phenotype.mechanics.maxAttachments = (int)parameter;
        else if( index == damage_rate_ind )
            cell.phenotype.cellInteractions.damageRate = parameter;
        else
            throw new IllegalArgumentException( "Incorrect signal with index " + index );
    }

    public void setSingleBehavior(Cell cell, String name, double parameter) throws Exception
    {
        int index = findBehaviorIndex( name );
        setSingleBehavior( cell, index, parameter );
    }

    public int findBehaviorIndex(String responseName) throws IllegalArgumentException
    {
        Integer result = behaviorToInt.get( responseName );
        if( result == null )
            throw new IllegalArgumentException( "Can not find behavior " + responseName );
        return result;
    }

    public double getSinglBehavior(Cell cell, String name) throws Exception
    {
        return getSingleBehavior( cell, findBehaviorIndex( name ) );
    }

    public double getSingleBehavior(Cell cell, int index) throws Exception
    {
        Microenvironment microenvironment = cell.getMicroenvironment();
        int m = microenvironment.numberDensities();
        int n = cell.getModel().getDefinitionsCount();

        if( index >= first_secretion_index && index < first_secretion_index + m )
            return cell.phenotype.secretion.secretionRates[index - first_secretion_index];
        if( index >= first_secretion_target_index && index < first_secretion_target_index + m )
            return cell.phenotype.secretion.saturationDensities[index - first_secretion_target_index];
        if( index >= first_uptake_index && index < first_uptake_index + m )
            return cell.phenotype.secretion.uptakeRates[index - first_uptake_index];
        if( index >= first_export_index && index < first_export_index + m )
            return cell.phenotype.secretion.netExportRates[index - first_export_index];

        int max_cycle_index = cell.phenotype.cycle.phases.size();
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
                return cell.phenotype.cycle.data.getExitRate( ind );
            return 0.0;
        }

        if( index == apoptosis_parameter_index )
            return cell.phenotype.death.rates.get( cell.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model ) );
        if( index == necrosis_parameter_index )
            return cell.phenotype.death.rates.get( cell.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model ) );
        if( index == migration_speed_index )
            return cell.phenotype.motility.migrationSpeed;
        if( index == migration_bias_index )
            return cell.phenotype.motility.migrationBias;
        if( index == persistence_time_index )
            return cell.phenotype.motility.persistenceTime;
        if( index >= first_chemotaxis_index && index < first_chemotaxis_index + m )
            return cell.phenotype.motility.chemotacticSensitivities[index - first_chemotaxis_index];
        if( index == cca_index )
            return cell.phenotype.mechanics.cellCellAdhesionStrength;
        if( index == elastic_index )
            return cell.phenotype.mechanics.attachmentElasticConstant;
        if( index >= first_affinity_index && index < first_affinity_index + n )
            return cell.phenotype.mechanics.cellAdhesionAffinities[index - first_affinity_index];
        if( index == max_adh_distance_index )
            return cell.phenotype.mechanics.relMaxAdhesionDistance;
        if( index == ccr_index )
            return cell.phenotype.mechanics.cellCellRepulsionStrength;
        if( index == cba_index )
            return cell.phenotype.mechanics.cellBMAdhesionStrength;
        if( index == cbr_index )
            return cell.phenotype.mechanics.cellBMRepulsionStrength;
        if( index == dead_phago_index )
            return cell.phenotype.cellInteractions.deadPhagocytosisRate;
        if( index >= first_phagocytosis_index && index < first_phagocytosis_index + n )
            return cell.phenotype.cellInteractions.livePhagocytosisRates[index - first_phagocytosis_index];
        if( index >= first_attack_index && index < first_attack_index + n )
            return cell.phenotype.cellInteractions.attackRates[index - first_attack_index];
        if( index >= first_fusion_index && index < first_fusion_index + n )
            return cell.phenotype.cellInteractions.fusionRates[index - first_fusion_index];
        if( index >= first_transformation_index && index < first_transformation_index + n )
            return cell.phenotype.cellTransformations.transformationRates[index - first_transformation_index];
        if( first_custom_ind >= 0 && index >= first_custom_ind && index < first_custom_ind + cell.customData.variables.size() )
            return cell.customData.variables.get( index - first_custom_ind ).value;
        if( index == movable_ind )
            return digitize( cell.isMovable );
        if( first_immunogenicity_index > -1 && index >= first_immunogenicity_index && index < first_immunogenicity_index + n )
            return cell.phenotype.cellInteractions.immunogenicities[index - first_immunogenicity_index];
        if( index == attachment_rate_ind )
            return cell.phenotype.mechanics.attachmentRate;
        if( index == detachment_rate_ind )
            return cell.phenotype.mechanics.detachmentRate;
        if( index == max_attachments_ind )
            return cell.phenotype.mechanics.maxAttachments;
        if( index == damage_rate_ind )
            return cell.phenotype.cellInteractions.damageRate;
        throw new Exception( "Warning: attempted to get behavior with unknown index!" );
    }

    public double[] getBaseBehaviors(Cell pCell)
    {
        CellDefinition cd = pCell.getModel().getCellDefinition( pCell.typeName );
        Model model = pCell.getModel();

        double[] parameters = new double[indexToBehavior.size()];

        insert( cd.phenotype.secretion.secretionRates, parameters, first_secretion_index );
        insert( cd.phenotype.secretion.saturationDensities, parameters, first_secretion_target_index );
        insert( cd.phenotype.secretion.uptakeRates, parameters, first_uptake_index );
        insert( cd.phenotype.secretion.netExportRates, parameters, first_export_index );

        int max_cycle_index = cd.phenotype.cycle.phases.size();
        if( max_cycle_index > 6 )
        {
            max_cycle_index = 6;
            System.out.println( "Warning: Standardized behaviors only support exit rate from the first 6 phases of a cell cycle!"
                    + "         Ignoring any later phase exit rates." );
        }
        for( int i = 0; i < max_cycle_index; i++ )
            parameters[first_cycle_index + i] = cd.phenotype.cycle.data.getExitRate( i );

        int apoptosis_index = cd.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );
        int necrosis_index = cd.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );

        parameters[apoptosis_parameter_index] = cd.phenotype.death.rates.get( apoptosis_index );
        parameters[necrosis_parameter_index] = cd.phenotype.death.rates.get( necrosis_index );
        parameters[migration_speed_index] = cd.phenotype.motility.migrationSpeed;
        parameters[migration_bias_index] = cd.phenotype.motility.migrationBias;
        parameters[persistence_time_index] = cd.phenotype.motility.persistenceTime;
        insert( cd.phenotype.motility.chemotacticSensitivities, parameters, first_chemotaxis_index );
        parameters[cca_index] = cd.phenotype.mechanics.cellCellAdhesionStrength;
        parameters[elastic_index] = cd.phenotype.mechanics.attachmentElasticConstant;
        String search_for1 = "adhesive affinity to " + model.getCellDefinition( 0 ).name;
        int first_affinity_index = findBehaviorIndex( search_for1 );
        insert( cd.phenotype.mechanics.cellAdhesionAffinities, parameters, first_affinity_index );
        parameters[max_adh_distance_index] = cd.phenotype.mechanics.relMaxAdhesionDistance;
        parameters[ccr_index] = cd.phenotype.mechanics.cellCellRepulsionStrength;
        parameters[cba_index] = cd.phenotype.mechanics.cellBMAdhesionStrength;
        parameters[cbr_index] = cd.phenotype.mechanics.cellBMRepulsionStrength;
        parameters[dead_phago_index] = cd.phenotype.cellInteractions.deadPhagocytosisRate;
        insert( cd.phenotype.cellInteractions.livePhagocytosisRates, parameters, first_phagocytosis_index );
        insert( cd.phenotype.cellInteractions.attackRates, parameters, first_attack_index );
        insert( cd.phenotype.cellInteractions.fusionRates, parameters, first_fusion_index );
        insert( cd.phenotype.cellTransformations.transformationRates, parameters, first_transformation_index );

        int first_custom_ind = findBehaviorIndex( "custom 0" );
        if( first_custom_ind >= 0 )
        {
            for( int nc = 0; nc < pCell.customData.variables.size(); nc++ )
                parameters[first_custom_ind + nc] = cd.custom_data.variables.get( nc ).value;
        }

        parameters[movable_ind] = digitize( cd.isMovable );
        insert( cd.phenotype.cellInteractions.immunogenicities, parameters, first_immunogenicity_index );
        parameters[attachment_rate_ind] = cd.phenotype.mechanics.attachmentRate;
        parameters[detachment_rate_ind] = cd.phenotype.mechanics.detachmentRate;
        parameters[max_attachments_ind] = cd.phenotype.mechanics.maxAttachments;
        parameters[damage_rate_ind] = cd.phenotype.cellInteractions.damageRate;
        return parameters;
    }

    private void insert(double[] source, double[] target, int pos)
    {
        System.arraycopy( source, 0, target, pos, source.length );
    }

    public double getSingleBaseBehavior(Cell cell, int index)
    {
        Microenvironment microenvironment = cell.getMicroenvironment();
        Model model = cell.getModel();
        int m = microenvironment.numberDensities();
        int n = model.getDefinitionsCount();

        CellDefinition cd = model.getCellDefinition( cell.typeName );

        if( index < 0 )


            if( index >= first_secretion_index && index < first_secretion_index + m )
                return cd.phenotype.secretion.secretionRates[index - first_secretion_index];

        if( index >= first_secretion_target_index && index < first_secretion_target_index + m )
            return cd.phenotype.secretion.saturationDensities[index - first_secretion_target_index];

        if( index >= first_uptake_index && index < first_uptake_index + m )
            return cd.phenotype.secretion.uptakeRates[index - first_uptake_index];

        if( index >= first_export_index && index < first_export_index + m )
            return cd.phenotype.secretion.netExportRates[index - first_export_index];

        int max_cycle_index = cd.phenotype.cycle.phases.size();
        if( max_cycle_index > 6 )
            throw new IllegalArgumentException(
                    "Warning: Standardized behaviors only support exit rate from the first 6 phases of a cell cycle!" );

        if( index >= first_cycle_index && index < first_cycle_index + 6 )
        {
            int ind = index - first_cycle_index;
            if( ind < max_cycle_index )
                return cd.phenotype.cycle.data.getExitRate( ind );
            return 0.0;
        }

        int apoptosis_index = cd.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );
        int necrosis_index = cd.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );

        if( index == apoptosis_parameter_index )
            return cd.phenotype.death.rates.get( apoptosis_index );
        if( index == necrosis_parameter_index )
            return cd.phenotype.death.rates.get( necrosis_index );
        if( index == migration_speed_index )
            return cd.phenotype.motility.migrationSpeed;
        if( index == migration_bias_index )
            return cd.phenotype.motility.migrationBias;
        if( index == persistence_time_index )
            return cd.phenotype.motility.persistenceTime;
        if( index >= first_chemotaxis_index && index < first_chemotaxis_index + m )
            return cd.phenotype.motility.chemotacticSensitivities[index - first_chemotaxis_index];
        if( index == cca_index )
            return cd.phenotype.mechanics.cellCellAdhesionStrength;
        if( index == elastic_index )
            return cd.phenotype.mechanics.attachmentElasticConstant;
        if( index >= first_affinity_index && index < first_affinity_index + n )
            return cd.phenotype.mechanics.cellAdhesionAffinities[index - first_affinity_index];
        if( index == max_adh_distance_index )
            return cd.phenotype.mechanics.relMaxAdhesionDistance;
        if( index == ccr_index )
            return cd.phenotype.mechanics.cellCellRepulsionStrength;
        if( index == cba_index )
            return cd.phenotype.mechanics.cellBMAdhesionStrength;
        if( index == cbr_index )
            return cd.phenotype.mechanics.cellBMRepulsionStrength;
        if( index == dead_phago_index )
            return cd.phenotype.cellInteractions.deadPhagocytosisRate;
        if( index >= first_phagocytosis_index && index < first_phagocytosis_index + n )
            return cd.phenotype.cellInteractions.livePhagocytosisRates[index - first_phagocytosis_index];
        if( index >= first_attack_index && index < first_attack_index + n )
            return cd.phenotype.cellInteractions.attackRates[index - first_attack_index];
        if( index >= first_fusion_index && index < first_fusion_index + n )
            return cd.phenotype.cellInteractions.fusionRates[index - first_fusion_index];
        if( index >= first_transformation_index && index < first_transformation_index + n )
            return cd.phenotype.cellTransformations.transformationRates[index - first_transformation_index];
        if( first_custom_ind >= 0 && index >= first_custom_ind && index < first_custom_ind + cell.customData.variables.size() )
            return cd.custom_data.variables.get( index - first_custom_ind ).value;
        if( index == movable_ind )
            digitize( cd.isMovable );
        if( first_immunogenicity_index > -1 && index >= first_immunogenicity_index && index < first_immunogenicity_index + n )
            return cd.phenotype.cellInteractions.immunogenicities[index - first_immunogenicity_index];
        if( index == attachment_rate_ind )
            return cd.phenotype.mechanics.attachmentRate;
        if( index == detachment_rate_ind )
            return cd.phenotype.mechanics.detachmentRate;
        if( index == max_attachments_ind )
            return cd.phenotype.mechanics.maxAttachments;
        if( index == damage_rate_ind )
            return cd.phenotype.cellInteractions.damageRate;
        throw new IllegalArgumentException( "Warning: attempted to get behavior with unknown index !" );
    }

    public double getSingleBaseBehavior(Model model, CellDefinition pCD, String name)
    {
        return getSingleBaseBehavior( model, pCD, findBehaviorIndex( name ) );
    }

    public double getSingleBaseBehavior(Model model, CellDefinition cd, int index)
    {
        Microenvironment microenvironment = cd.getMicroenvironment();
        int m = microenvironment.numberDensities();
        int n = model.getDefinitionsCount();

        if( index >= first_secretion_index && index < first_secretion_index + m )
            return cd.phenotype.secretion.secretionRates[index - first_secretion_index];
        if( index >= first_secretion_target_index && index < first_secretion_target_index + m )
            return cd.phenotype.secretion.saturationDensities[index - first_secretion_target_index];
        if( index >= first_uptake_index && index < first_uptake_index + m )
            return cd.phenotype.secretion.uptakeRates[index - first_uptake_index];
        if( index >= first_export_index && index < first_export_index + m )
            return cd.phenotype.secretion.netExportRates[index - first_export_index];

        int max_cycle_index = cd.phenotype.cycle.phases.size();
        if( max_cycle_index > 6 )
            throw new IllegalArgumentException(
                    "Warning: Standardized behaviors only support exit rate from the first 6 phases of a cell cycle!" );
        if( index >= first_cycle_index && index < first_cycle_index + 6 )
        {
            int ind = index - first_cycle_index;
            if( ind < max_cycle_index )
                return cd.phenotype.cycle.data.getExitRate( ind );
            return 0.0;
        }

        int apoptosis_index = cd.phenotype.death.findDeathModelIndex( PhysiCellConstants.apoptosis_death_model );
        int necrosis_index = cd.phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );

        if( index == apoptosis_parameter_index )
            return cd.phenotype.death.rates.get( apoptosis_index );
        if( index == necrosis_parameter_index )
            return cd.phenotype.death.rates.get( necrosis_index );
        if( index == migration_speed_index )
            return cd.phenotype.motility.migrationSpeed;
        if( index == migration_bias_index )
            return cd.phenotype.motility.migrationBias;
        if( index == persistence_time_index )
            return cd.phenotype.motility.persistenceTime;
        if( index >= first_chemotaxis_index && index < first_chemotaxis_index + m )
            return cd.phenotype.motility.chemotacticSensitivities[index - first_chemotaxis_index];
        if( index == cca_index )
            return cd.phenotype.mechanics.cellCellAdhesionStrength;
        if( index == elastic_index )
            return cd.phenotype.mechanics.attachmentElasticConstant;
        if( index >= first_affinity_index && index < first_affinity_index + n )
            return cd.phenotype.mechanics.cellAdhesionAffinities[index - first_affinity_index];
        if( index == max_adh_distance_index )
            return cd.phenotype.mechanics.relMaxAdhesionDistance;
        if( index == ccr_index )
            return cd.phenotype.mechanics.cellCellRepulsionStrength;
        if( index == cba_index )
            return cd.phenotype.mechanics.cellBMAdhesionStrength;
        if( index == cbr_index )
            return cd.phenotype.mechanics.cellBMRepulsionStrength;
        if( index == dead_phago_index )
            return cd.phenotype.cellInteractions.deadPhagocytosisRate;
        if( index >= first_phagocytosis_index && index < first_phagocytosis_index + n )
            return cd.phenotype.cellInteractions.livePhagocytosisRates[index - first_phagocytosis_index];
        if( index >= first_attack_index && index < first_attack_index + n )
            return cd.phenotype.cellInteractions.attackRates[index - first_attack_index];
        if( index >= first_fusion_index && index < first_fusion_index + n )
            return cd.phenotype.cellInteractions.fusionRates[index - first_fusion_index];
        if( index >= first_transformation_index && index < first_transformation_index + n )
            return cd.phenotype.cellTransformations.transformationRates[index - first_transformation_index];
        if( first_custom_ind >= 0 && index >= first_custom_ind && index < first_custom_ind + cd.custom_data.variables.size() )
            return cd.custom_data.variables.get( index - first_custom_ind ).value;
        if( index == movable_ind )
            return digitize( cd.isMovable );
        if( first_immunogenicity_index > -1 && index >= first_immunogenicity_index && index < first_immunogenicity_index + n )
            return cd.phenotype.cellInteractions.immunogenicities[index - first_immunogenicity_index];
        if( index == attachment_rate_ind )
            return cd.phenotype.mechanics.attachmentRate;
        if( index == detachment_rate_ind )
            return cd.phenotype.mechanics.detachmentRate;
        if( index == max_attachments_ind )
            return cd.phenotype.mechanics.maxAttachments;
        if( index == damage_rate_ind )
            return cd.phenotype.cellInteractions.damageRate;
        throw new IllegalArgumentException( "Warning: attempted to get behavior with unknown index !" );
    }

    public double getSinglBaseBehavior(Cell cell, String name)
    {
        return getSingleBaseBehavior( cell, findBehaviorIndex( name ) );
    }

    public double[] getBaseBehaviors(Cell cell, int[] indices)
    {
        double[] parameters = new double[indices.length];
        for( int n = 0; n < indices.length; n++ )
            parameters[n] = getSingleBaseBehavior( cell, indices[n] );
        return parameters;
    }

    public double[] getBaseBehaviors(Cell cell, String[] names)
    {
        double[] parameters = new double[names.length];
        for( int n = 0; n < names.length; n++ )
            parameters[n] = getSinglBaseBehavior( cell, names[n] );
        return parameters;
    }
}