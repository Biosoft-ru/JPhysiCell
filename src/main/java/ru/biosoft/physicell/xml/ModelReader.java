package ru.biosoft.physicell.xml;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.MicroenvironmentOptions;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.CellCSVReader;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellInteractions;
import ru.biosoft.physicell.core.CellTransformations;
import ru.biosoft.physicell.core.CycleModel;
import ru.biosoft.physicell.core.Death;
import ru.biosoft.physicell.core.DeathParameters;
import ru.biosoft.physicell.core.Mechanics;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Motility;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.PhysiCellSettings;
import ru.biosoft.physicell.core.Secretion;
import ru.biosoft.physicell.core.StandardModels;
import ru.biosoft.physicell.core.Volume;

public class ModelReader extends Constants
{

    public Model read(File f) throws Exception
    {
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document doc = builder.parse( f );
        NodeList nodes = doc.getChildNodes();

        Element physicell = findElement( nodes, PHYSICELL_ELEMENT );
        if( physicell == null )
            throw new Exception( "Physicell base element not found" );

        Model model = new Model();
        Microenvironment m = model.getMicroenvironment();
        readDomain( physicell, m );
        readOverall( physicell, model );
        readOptions( physicell );
        readSave( physicell, model );
        readMicroenvironmentSetup( physicell, m );
        readCellDefinitions( physicell, m );
        readPhenotypes( physicell, m );
        readInitialConditions( physicell, model );
        readUserParameters( physicell, model );

        return model;
    }

    private void readDomain(Element physicell, Microenvironment m)
    {
        Element domainElement = findElement( physicell, DOMAIN_ELEMENT );
        double minX = 0;
        double maxX = 100;
        double minY = 0;
        double maxY = 100;
        double minZ = 0;
        double maxZ = 100;
        double dx = 10;
        double dy = 10;
        double dz = 10;
        boolean use2D = false;
        for( Element el : getAllElements( domainElement ) )
        {
            switch( el.getTagName() )
            {
                case X_MIN:
                    minX = getDoubleVal( el );
                    break;
                case X_MAX:
                    maxX = getDoubleVal( el );
                    break;
                case Y_MIN:
                    minY = getDoubleVal( el );
                    break;
                case Y_MAX:
                    maxY = getDoubleVal( el );
                    break;
                case Z_MIN:
                    minZ = getDoubleVal( el );
                    break;
                case Z_MAX:
                    maxZ = getDoubleVal( el );
                    break;
                case DX:
                    dx = getDoubleVal( el );
                    break;
                case DY:
                    dy = getDoubleVal( el );
                    break;
                case DZ:
                    dz = getDoubleVal( el );
                    break;
                case USE_2D:
                    use2D = getBoolVal( el );
                    break;
            }
        }
        m.resizeSpace( minX, maxX, minY, maxY, minZ, maxZ, dx, dy, dz );
        //        efault_microenvironment_options.simulate_2D = xml_get_bool_value( node, "use_2D" ); 
        //
        if( m.options.simulate_2D )
        {
            minZ = -0.5 * dz;
            maxZ = 0.5 * dz;
        }
        m.options.X_range = new double[] {minX, maxX};
        m.options.Y_range = new double[] {minY, maxY};
        m.options.Z_range = new double[] {minZ, maxZ};

        m.options.dx = dx;
        m.options.dy = dy;
        m.options.dz = dz;
        m.options.simulate_2D = use2D;
    }

    public void readOverall(Element physicell, Model model)
    {
        Microenvironment m = model.getMicroenvironment();
        double maxTime = 100;
        //        String maxTimeUnits = "min";
        String timeUnits = "min";
        String spaceUnits = "micron";
        //        String diffusionUnits = "min";
        //        String mechanicsUnits = "min";
        //        String phenotypeUnits = "min";
        double diffusionStep = 0.01;
        double mechanicsStep = 0.1;
        double phenotypeSteps = 6;
        Element overallElement = findElement( physicell, OVERALL );
        for( Element el : getAllElements( overallElement ) )
        {
            switch( el.getTagName() )
            {
                case MAX_TIME:
                    maxTime = getDoubleVal( el );
                    //                    maxTimeUnits = el.getAttribute( UNITS );
                    break;
                case TIME_UNITS:
                    timeUnits = getVal( el );
                    break;
                case SPACE_UNITS:
                    spaceUnits = getVal( el );
                    break;
                case DT_DIFFUSION:
                    //                    mechanicsUnits = el.getAttribute( UNITS );
                    diffusionStep = getDoubleVal( el );
                    break;
                case DT_MECHANICS:
                    //                    diffusionUnits = el.getAttribute( UNITS );
                    mechanicsStep = getDoubleVal( el );
                    break;
                case DT_PHENOTYPE:
                    //                    phenotypeUnits = el.getAttribute( UNITS );
                    phenotypeSteps = getDoubleVal( el );
                    break;

            }
        }
        m.timeUnits = timeUnits;
        m.spatialUnits = spaceUnits;
        model.setDiffusionDt( diffusionStep );
        model.setPhenotypeDt( phenotypeSteps );
        model.setMechanicsDt( mechanicsStep );
        model.setTMax( maxTime );
    }


    public void readSave(Element physicell, Model model)
    {
        Element saveElement = findElement( physicell, "save" );
        for( Element el : getAllElements( saveElement ) )
        {
            switch( el.getTagName() )
            {
                case "folder":
                    String folderName = getVal( el );
                    break;
                case "full_data":
                    Element intervalElement = findElement( el, "interval" );
                    Element enabledElement = findElement( el, "enable" );
                    if( enabledElement != null )
                    {
                        model.setEnableFullSaves( getBoolVal( enabledElement ) );
                    }
                    if( intervalElement != null && model.isEnableFullSaves() )
                    {
                        double interval = getDoubleVal( intervalElement );
                        String intervalUnits = getAttr( intervalElement, "units" );
                        model.setSaveInterval( interval );
                    }

                    break;
                case "SVG":
                    intervalElement = findElement( el, "interval" );
                    enabledElement = findElement( el, "enable" );
                    boolean enabled = false;
                    if( enabledElement != null )
                    {
                        enabled = getBoolVal( enabledElement );
                    }
                    if( intervalElement != null && enabled )
                    {
                        double interval = getDoubleVal( intervalElement );
                        String intervalUnits = getAttr( intervalElement, "units" );
                    }

                    break;
                case "legacy_data":
                    enabledElement = findElement( el, "" );
                    if( enabledElement != null )
                    {
                        enabled = getBoolVal( enabledElement );
                    }
                    break;
            }
        }
    }

    public void readOptions(Element physicell)
    {
        Element optionsElement = findElement( physicell, "options" );
        for( Element el : getAllElements( optionsElement ) )
        {
            switch( el.getTagName() )
            {
                //                case "legacy_random_points_on_sphere_in_divide":
                //                    boolean legacy_random_points_on_sphere_in_divide = getBoolVal( el );
                //                    PhysiCellSettings.cell_division_orientation = LegacyRandomOnUnitSphere;
                //                    break;
                case "virtual_wall_at_domain_edge":
                {
                    try
                    {
                        boolean virtual_wall_at_domain_edge = getBoolVal( el );
                        if( virtual_wall_at_domain_edge )
                            StandardModels
                                    .getDefaultCellDefinition().functions.add_cell_basement_membrane_interactions = new StandardModels.standard_domain_edge_avoidance_interactions();
                        break;
                    }
                    catch( Exception ex )
                    {

                    }
                }
                case "disable_automated_spring_adhesions":
                    PhysiCellSettings.disable_automated_spring_adhesions = getBoolVal( el );
                    break;
            }
        }
    }

    public void readMicroenvironmentSetup(Element physicell, Microenvironment m)
    {
        Element microenvironmentSetupElement = findElement( physicell, "microenvironment_setup" );
        List<Element> variableElements = findAllElements( microenvironmentSetupElement, "variable" );

        boolean activated_Dirichlet_boundary_detected = false;

        int varSize = variableElements.size();
        double[] initialValues = new double[varSize];
        double[] Dirichlet_condition = new double[varSize];
        boolean[] Dirichlet_activation = new boolean[varSize];

        boolean[] Dirichlet_xmin = new boolean[varSize];
        boolean[] Dirichlet_xmax = new boolean[varSize];
        boolean[] Dirichlet_ymin = new boolean[varSize];
        boolean[] Dirichlet_ymax = new boolean[varSize];
        boolean[] Dirichlet_zmin = new boolean[varSize];
        boolean[] Dirichlet_zmax = new boolean[varSize];

        double[] Dirichlet_xmin_values = new double[varSize];
        double[] Dirichlet_xmax_values = new double[varSize];
        double[] Dirichlet_ymin_values = new double[varSize];
        double[] Dirichlet_ymax_values = new double[varSize];
        double[] Dirichlet_zmin_values = new double[varSize];
        double[] Dirichlet_zmax_values = new double[varSize];
        boolean[] Dirichlet_all = new boolean[varSize];

        for( int i = 0; i < varSize; i++ )
        {
            Element variableElement = variableElements.get( i );
            String name = getAttr( variableElement, "name" );
            String units = getAttr( variableElement, "units" );
            int id = getIntAttr( variableElement, "ID" );
            Element physicalElement = findElement( variableElement, "physical_parameter_set" );
            Element diffusionElement = findElement( physicalElement, "diffusion_coefficient" );
            //            String diffusionUnits = getAttr( diffusionElement, "units" );
            double diffusionValue = getDoubleVal( diffusionElement );
            Element decayElement = findElement( physicalElement, "decay_rate" );
            //            String decayUnits = getAttr( diffusionElement, "units" );
            double decayValue = getDoubleVal( decayElement );

            if( id == 0 )
                m.setDensity( id, name, units, diffusionValue, decayValue );
            else
                m.addDensity( name, units, diffusionValue, decayValue );

            Element initialConditionElement = findElement( variableElement, "initial_condition" );
            //            String initialUnits = initialConditionElement.getAttribute( "units" );
            initialValues[i] = getDoubleVal( initialConditionElement );

            Element DirichletBoundaryElement = findElement( variableElement, "Dirichlet_boundary_condition" );
            //            String dirichletUnits = getAttr( DirichletBoundaryElement, "units" );
            boolean dirichletEnabled = getBoolAttr( DirichletBoundaryElement, "enabled" );
            double dirichletValue = getDoubleVal( DirichletBoundaryElement );
            Dirichlet_condition[i] = dirichletValue;
            Dirichlet_activation[i] = dirichletEnabled;
            Dirichlet_all[i] = dirichletEnabled;

            if( dirichletEnabled )
                activated_Dirichlet_boundary_detected = true;

            Dirichlet_xmin[i] = dirichletEnabled;
            Dirichlet_xmax[i] = dirichletEnabled;
            Dirichlet_ymin[i] = dirichletEnabled;
            Dirichlet_ymax[i] = dirichletEnabled;
            Dirichlet_zmin[i] = dirichletEnabled;
            Dirichlet_zmax[i] = dirichletEnabled;

            Dirichlet_xmin_values[i] = dirichletValue;
            Dirichlet_xmax_values[i] = dirichletValue;
            Dirichlet_ymin_values[i] = dirichletValue;
            Dirichlet_ymax_values[i] = dirichletValue;
            Dirichlet_zmin_values[i] = dirichletValue;
            Dirichlet_zmax_values[i] = dirichletValue;

            Element dirichletOptionsElement = findElement( variableElement, "Dirichlet_options" );
            if( dirichletOptionsElement != null )
            {
                for( Element el : findAllElements( dirichletOptionsElement, "boundary_value" ) )
                {
                    boolean enabled = getBoolAttr( el, "enabled" );
                    double value = getDoubleVal( el );
                    if( !enabled )
                        Dirichlet_all[i] = false;
                    switch( getAttr( el, "ID" ) )
                    {
                        case "xmin":
                            Dirichlet_xmin[i] = enabled;
                            Dirichlet_xmin_values[i] = value;
                            break;
                        case "xmax":
                            Dirichlet_xmax[i] = enabled;
                            Dirichlet_xmax_values[i] = value;
                            break;
                        case "ymin":
                            Dirichlet_ymin[i] = enabled;
                            Dirichlet_ymin_values[i] = value;
                            break;
                        case "ymax":
                            Dirichlet_ymax[i] = enabled;
                            Dirichlet_ymax_values[i] = value;
                            break;
                        case "zmin":
                            Dirichlet_zmin[i] = enabled;
                            Dirichlet_zmin_values[i] = value;
                            break;
                        case "zmax":
                            Dirichlet_zmax[i] = enabled;
                            Dirichlet_zmax_values[i] = value;
                            break;
                    }
                }
            }
        }
        MicroenvironmentOptions options = m.getOptions();
        options.outer_Dirichlet_conditions = false;

        options.Dirichlet_condition_vector = Dirichlet_condition;
        options.Dirichlet_activation_vector = Dirichlet_activation;
        options.initial_condition_vector = initialValues;

        options.Dirichlet_all = Dirichlet_all;

        options.Dirichlet_xmin = Dirichlet_xmin;
        options.Dirichlet_xmax = Dirichlet_xmax;
        options.Dirichlet_ymin = Dirichlet_ymin;
        options.Dirichlet_ymax = Dirichlet_ymax;
        options.Dirichlet_zmin = Dirichlet_zmin;
        options.Dirichlet_zmax = Dirichlet_zmax;

        options.Dirichlet_xmin_values = Dirichlet_xmin_values;
        options.Dirichlet_xmax_values = Dirichlet_xmax_values;
        options.Dirichlet_ymin_values = Dirichlet_ymin_values;
        options.Dirichlet_ymax_values = Dirichlet_ymax_values;
        options.Dirichlet_zmin_values = Dirichlet_zmin_values;
        options.Dirichlet_zmax_values = Dirichlet_zmax_values;

        for( int i = 0; i < m.number_of_voxels(); i++ )
            m.density[i] = options.initial_condition_vector.clone();

        // if any of the substrates have outer Dirichlet conditions enables, then set the outer_Dirichlet_conditions = true;       
        if( activated_Dirichlet_boundary_detected )
            options.outer_Dirichlet_conditions = true;

        Element calculate_gradientsElement = findElement( microenvironmentSetupElement, "calculate_gradients" );
        if( calculate_gradientsElement != null )
            options.calculate_gradients = getBoolVal( calculate_gradientsElement );

        Element trackSubstrate = findElement( microenvironmentSetupElement, "track_internalized_substrates_in_each_agent" );
        if( trackSubstrate != null )
            options.track_internalized_substrates_in_each_agent = getBoolVal( trackSubstrate );

        Microenvironment.initialize_microenvironment( m );
    }

    //    void initialize_cell_definitions_from_pugixml(Element el)
    //    {
    //        //        pugi::xml_node node_options; 
    //
    //        //        node_options = xml_find_node( root , "options" ); 
    //        //        if( node_options )
    //        //        {
    //        //            bool settings = 
    //        //                xml_get_bool_value( node_options, "virtual_wall_at_domain_edge" ); 
    //        //            if( settings )
    //        //            {
    //        //                std::cout << "virtual_wall_at_domain_edge: enabled" << std::endl; 
    //        //                cell_defaults.functions.add_cell_basement_membrane_interactions = standard_domain_edge_avoidance_interactions;
    //        //            }
    //        //
    //        //        }
    //        //        
    //        //        // first, let's pre-build the map. 
    //        //        prebuild_cell_definition_index_maps(); 
    //        //        
    //        Element cdsElement = findElement( el, "cell_definitions" );
    //
    //        for( Element cdElement : findAllElements( cdsElement, "cell_definition" ) )
    //        {
    //            System.out.println( "Processing " + getAttr( cdElement, "name" ) + " ... " );
    //
    //            initialize_cell_definition_from_pugixml( cdElement );
    //            build_cell_definitions_maps();
    //        }
    //    }
    //    Map<String, Integer> cellDefinitions = null;
    //    int cellDefinitionSize = 0;
    //    public void preReadCellDefinitions(Element physicell, Microenvironment m) throws Exception
    //    {
    //        //        cellDefinitions = new HashMap<>();
    //        Element cellDefinitionsElement = findElement( physicell, "cell_definitions" );
    //        for( Element cdElement : findAllElements( cellDefinitionsElement, "cell_definition" ) )
    //        {
    //            cellDefinitionSize++;
    //            //            String name = getAttr( cdElement, "name" );
    //            //            Integer ID = getIntAttr( cdElement, "ID" );
    //            //            cellDefinitions.put( name, ID );
    //        }
    //    }

    public void readCellDefinitions(Element physicell, Microenvironment m) throws Exception
    {
        Element cellDefinitionsElement = findElement( physicell, "cell_definitions" );
        for( Element cdElement : findAllElements( cellDefinitionsElement, "cell_definition" ) )
        {
            CellDefinition cd;
            String name = getAttr( cdElement, "name" );
            Integer ID = getIntAttr( cdElement, "ID" );
            int id = ID == null ? -1 : ID;
            boolean defaultDefinition = name.equals( "default" );// || id == 0;
            if( defaultDefinition ) //TODO: check ID=0 
                cd = StandardModels.createFromDefault( name, id, m );
            else
                cd = new CellDefinition( m, id, name );

            CellDefinition parent = null;
            String parentType = getAttr( cdElement, "parent_type" );
            //                        if( ! )
            //                            parent = CellDefinition.getCellDefinition( parentType );
            //            boolean use_default_as_parent_without_specifying = false;
            if( parentType.isEmpty() && !defaultDefinition )
            {
                parent = StandardModels.createFromDefault( "default", -1, m );
                cd = parent.clone( name, id, m );
                //                use_default_as_parent_without_specifying = true;
            }
            //
            //            // if we found something to inherit from, then do it! 
            //            if( parent != null )
            //            {
            //                System.out.println( "\tCopying from type " + parent.name + " ... " );
            //                cd = parent.clone( name, id, m );
            //            }
            CellDefinition.registerCellDefinition( cd );
        }
    }

    //    public void readCellDefinitions2(Element physicell, Microenvironment m) throws Exception
    //    {
    //        Element cellDefinitionsElement = findElement( physicell, "cell_definitions" );
    //        for( Element cdElement : findAllElements( cellDefinitionsElement, "cell_definition" ) )
    //        {
    //            CellDefinition cd;
    //            String name = getAttr( cdElement, "name" );
    //            Integer ID = getIntAttr( cdElement, "ID" );
    //            int id = ID == null ? -1 : ID;
    //            boolean defaultDefinition = name.equals( "default" ) || id == 0;
    //            if( defaultDefinition ) //TODO: check ID=0 
    //                cd = StandardModels.createDefaultCellDefinition( name, m );
    //            else
    //                cd = new CellDefinition( m, name );
    //            cd.type = id;
    //            CellDefinition.registerCellDefinition( cd );
    //
    //            CellDefinition parent = null;
    //            String parentType = getAttr( cdElement, "parent_type" );
    //            if( !parentType.isEmpty() )
    //                parent = CellDefinition.getCellDefinition( parentType );
    //            //            boolean use_default_as_parent_without_specifying = false;
    //            if( parent == null && !defaultDefinition )
    //            {
    //                parent = StandardModels.createDefaultCellDefinition( "default", m );
    //                //                use_default_as_parent_without_specifying = true;
    //            }
    //
    //            // if we found something to inherit from, then do it! 
    //            if( parent != null )
    //            {
    //                System.out.println( "\tCopying from type " + parent.name + " ... " );
    //                cd = parent;
    //                // but recover the name and ID (type)
    //                cd.name = name;
    //                cd.type = id;
    //            }
    public void readPhenotypes(Element physicell, Microenvironment m) throws Exception
    {


        Element cellDefinitionsElement = findElement( physicell, "cell_definitions" );
        for( Element cdElement : findAllElements( cellDefinitionsElement, "cell_definition" ) )
        {
            String name = getAttr( cdElement, "name" );
            CellDefinition cd = CellDefinition.getCellDefinition( name );
            Phenotype p = cd.phenotype;

            CellDefinition parent = null;
            String parentType = getAttr( cdElement, "parent_type" );
            if( !parentType.isEmpty() )
                parent = CellDefinition.getCellDefinition( parentType );
            //            boolean use_default_as_parent_without_specifying = false;
            //            if( parent == null && !defaultDefinition )
            //            {
            //                parent = StandardModels.createFromDefault( "default", -1, m );
            //                //                use_default_as_parent_without_specifying = true;
            //            }

            // if we found something to inherit from, then do it! 
            if( parent != null )
            {
                System.out.println( "\tCopying from type " + parent.name + " ... " );
                cd = parent.clone( name, cd.type, m );
                CellDefinition.registerCellDefinition( cd );
            }

            p.sync( m );
            Element phenotypeElement = findElement( cdElement, "phenotype" );
            for( Element el : getAllElements( phenotypeElement ) )
            {
                switch( el.getTagName() )
                {
                    case "cycle":
                        readCycle( el, cd );
                        break;
                    case "death":
                        readDeath( el, p );
                        break;
                    case "volume":
                        readVolume( el, p );
                        break;
                    case "mechanics":
                        readMechanics( el, p );
                        break;
                    case "motility":
                        readMotility( el, cd, m );
                        break;
                    case "secretion":
                        readSecretion( el, p, m );
                        break;
                    case "cell_interactions":
                        readCellInteractions( el, cd );
                        break;
                    case "cell_transformations":
                        readCellTransformations( el, cd );
                        break;
                }
            }
            Element customDataEl = findElement( cdElement, "custom_data" );
            readCustomData( customDataEl, cd );
        }
    }

    private void readCycle(Element el, CellDefinition cd) throws Exception
    {
        Phenotype p = cd.phenotype;
        CycleModel model = p.cycle;
        String name = getAttr( el, "name" );
        Integer code = getIntAttr( el, "code" );
        if( code == null )
            code = PhysiCellConstants.find_cycle_model_code( name );
        if( code < 0 )
            throw new Exception( "Error. Unable to identify cycle model\n" + "node.attribute(\"name\").value()\n" + "("
                    + getIntAttr( el, "code" ) + ")" );
        model.add_phase( code, name );

        // Set the model, but only if it was specified. 
        if( getIntAttr( el, "code" ) != null )
        {
            // set the model 
            //switch( model )   // do not use a switch stmt to avoid compile errors related to "static const int" on various compilers
            if( code == PhysiCellConstants.advanced_Ki67_cycle_model )
            {
                p.cycle = StandardModels.Ki67_advanced.clone();
            }
            else if( code == PhysiCellConstants.basic_Ki67_cycle_model )
            {
                p.cycle = StandardModels.Ki67_basic.clone();
            }
            else if( code == PhysiCellConstants.flow_cytometry_cycle_model )
            {
                p.cycle = StandardModels.flow_cytometry_cycle_model.clone();
            }
            else if( code == PhysiCellConstants.live_apoptotic_cycle_model ) // ?
            {
                p.cycle = StandardModels.live.clone(); // ?
                System.out.println( "Warning: live_apoptotic_cycle_model not directly supported.\n"
                        + "         Substituting live cells model. Set death rates=0." );
            }
            else if( code == PhysiCellConstants.total_cells_cycle_model )
            {
                p.cycle = StandardModels.live.clone();
                System.out.println( "Warning: total_cells_cycle_model not directly supported.\n"
                        + "         Substituting live cells model. Set death rates=0." );
            }
            else if( code == PhysiCellConstants.live_cells_cycle_model )
            {
                p.cycle = StandardModels.live.clone();
            }
            else if( code == PhysiCellConstants.flow_cytometry_separated_cycle_model )
            {
                p.cycle = StandardModels.flow_cytometry_separated_cycle_model.clone();
            }
            else if( code == PhysiCellConstants.cycling_quiescent_model )
            {
                p.cycle = StandardModels.cycling_quiescent.clone();
            }
            else
            {
                throw new Exception( "Warning: Unknown cycle model " );
            }
            //            cd.phenotype.cycle.sync_to_cycle_model( cd.functions.cycle_model );
        }

        Element transitionElement = findElement( el, "phase_transition_rates" );
        //        String units = getAttr( transitionElement, "units" );
        if( transitionElement != null )
        {
            for( Element rate : findAllElements( transitionElement, "rate" ) )
            {
                try
                {
                    int startIndex = getIntAttr( rate, "start_index" );
                    int endIndex = getIntAttr( rate, "end_index" );
                    boolean fixedDuration = getBoolAttr( rate, "fixed_duration" );
                    double rateValue = getDoubleVal( rate );
                    p.cycle.setTransitionRate( startIndex, endIndex, rateValue );
                    p.cycle.phase_link( startIndex, endIndex ).fixedDuration = fixedDuration;
                }
                catch( Exception ex )
                {
                    ex.printStackTrace();//TODO: change to something nice
                }
            }
        }
        // Check for phase durations (as an alternative to transition rates)
        Element durationElement = findElement( el, "phase_durations" );
        //        if( node.child( "phase_durations" ) )
        //        { node = node.child( "phase_durations" ); }
        if( durationElement != null )
        {
            for( Element duration : findAllElements( durationElement, "duration" ) )
            {
                int start = getIntAttr( duration, "index" );
                boolean fixed = getBoolAttr( duration, "fixed_duration" );
                // actual value of the duration 
                double value = getDoubleVal( duration );
                // set the transition rate 
                p.cycle.data.setExitRate( start, 1.0 / ( value + 1e-16 ) );
                // set it to fixed / non-fixed 
                p.cycle.phase_links.get( start ).get( 0 ).fixedDuration = fixed;
            }
        }

        //    
        // now, if we inherited from another cell, AND 
        // if that parent type has the same cylce model, 
        //  then overwrite with their transition rates 
        //    
        //    if( pParent != null )
        //    {
        //        if( cd.phenotype.cycle.code == pParent.phenotype.cycle.code )
        //        {
        //            System.out.println( "copying data ..." );
        //            System.out.println( pParent.name + " to " + cd.name );
        //            cd.phenotype.cycle.data = pParent.phenotype.cycle.data;
        //        }
        //    }
    }

    private void readDeath(Element el, Phenotype p) throws Exception
    {
        Death death = p.death;
        for( Element modelElement : findAllElements( el, "model" ) )
        {
            Integer code = getIntAttr( modelElement, "code" );
            String name = getAttr( modelElement, "name" );

            if( code == null )
                code = PhysiCellConstants.find_cycle_model_code( name );
            if( code < 0 )
                throw new Exception( "Can not identify death model " + name );

            int death_index = death.findDeathModelIndex( code );

            boolean alreadyExists = false;
            if( death.models.get( death_index ).code == code )
            {
                alreadyExists = true;
            }

            //            CycleModel model = new CycleModel();
            //            model.name = name;
            //            model.code = code;

            DeathParameters parameters = new DeathParameters();
            Element deathRateElement = findElement( modelElement, "death_rate" );
            double rateValue = getDoubleVal( deathRateElement );
            parameters.time_units = getAttr( deathRateElement, "units" );
            Element durationsElement = findElement( modelElement, "phase_durations" );
            if( durationsElement != null )
            {
                //            String durationUnits = getAttr( durationsElement, "units" );
                for( Element durationElement : findAllElements( durationsElement, "duration" ) )
                {
                    int index = getIntAttr( durationElement, "index" );
                    boolean fixedDuration = getBoolAttr( durationElement, "fixed_duration" );
                    double durationValue = getDoubleVal( durationElement );
                    p.death.models.get( death_index ).data.setExitRate( index, 1.0 / ( durationValue + 1e-16 ) );
                    p.death.models.get( death_index ).phase_links.get( index ).get( 0 ).fixedDuration = fixedDuration;
                }
            }

            Element transitionsElement = findElement( modelElement, "phase_transition_rates" );
            if( transitionsElement != null )
            {
                //            String durationUnits = getAttr( durationsElement, "units" );
                for( Element transitionElement : findAllElements( transitionsElement, "rate" ) )
                {
                    int start = getIntAttr( transitionElement, "start_index" );
                    int end = getIntAttr( transitionElement, "end_index" );
                    boolean fixed = getBoolAttr( transitionElement, "fixed_duration" );
                    double value = getDoubleVal( transitionElement );
                    p.death.models.get( death_index ).setTransitionRate( start, end, value );
                    p.death.models.get( death_index ).phase_link( start, end ).fixedDuration = fixed;
                }
            }

            Element parametersElement = findElement( modelElement, "parameters" );
            for( Element paramEl : getAllElements( parametersElement ) )
            {
                switch( paramEl.getTagName() )
                {
                    case "unlysed_fluid_change_rate":
                        parameters.unlysed_fluid_change_rate = getDoubleVal( paramEl );
                        //                        String unlysed_fluid_change_rate_units = getAttr( paramEl, "units" ); //TODO
                        break;
                    case "lysed_fluid_change_rate":
                        parameters.lysed_fluid_change_rate = getDoubleVal( paramEl );
                        //                        String lysed_fluid_change_rate_units = getAttr( paramEl, "units" );
                        break;
                    case "cytoplasmic_biomass_change_rate":
                        parameters.cytoplasmic_biomass_change_rate = getDoubleVal( paramEl );
                        //                        String cytoplasmic_biomass_change_rate_units = getAttr( paramEl, "units" );
                        break;
                    case "nuclear_biomass_change_rate":
                        parameters.nuclear_biomass_change_rate = getDoubleVal( paramEl );
                        //                        String nuclear_biomass_change_rate_units = getAttr( paramEl, "units" );
                        break;
                    case "calcification_rate":
                        parameters.calcification_rate = getDoubleVal( paramEl );
                        //                        String calcification_rate_units = getAttr( paramEl, "units" );
                        break;
                    case "relative_rupture_volume":
                        parameters.relative_rupture_volume = getDoubleVal( paramEl );
                        //                        String relative_rupture_volume_units = getAttr( paramEl, "units" );
                        break;
                }
            }

            if( code == PhysiCellConstants.apoptosis_death_model )
            {
                if( !alreadyExists )
                {
                    death.addDeathModel( rateValue, StandardModels.apoptosis, parameters );
                    death_index = death.findDeathModelIndex( code );
                }
                else
                {
                    death.parameters.set( death_index, parameters );
                    death.rates.set( death_index, rateValue );
                }
            }
            else if( code == PhysiCellConstants.necrosis_death_model ) // set necrosis parameters 
            {
                if( !alreadyExists )
                {
                    death.addDeathModel( rateValue, StandardModels.necrosis, parameters );
                    death_index = death.findDeathModelIndex( code );
                }
                else
                {
                    death.parameters.set( death_index, parameters );
                    death.rates.set( death_index, rateValue );
                }
            }
            else if( code == PhysiCellConstants.autophagy_death_model )
            {
                System.out.println( "Warning: autophagy_death_model not yet supported.<br>         Skipping this model." );
            }
            else
            {
                throw new Exception( "Warning: Unknown death model " );
            }
        }
        //        <death>
        //                <model code="100" name="apoptosis">
        //                  <death_rate units="1/min">0</death_rate>
        //                  <phase_durations units="min">
        //                    <duration index="0" fixed_duration="true">516</duration>
        //                  </phase_durations>
        //                  <parameters>
        //                    <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
        //                    <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
        //                    <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
        //                    <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
        //                    <calcification_rate units="1/min">0</calcification_rate>
        //                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
        //                  </parameters>
        //                </model>
        //                <model code="101" name="necrosis">
        //                  <death_rate units="1/min">0.0</death_rate>
        //                  <phase_durations units="min">
        //                    <duration index="0" fixed_duration="true">0</duration>
        //                  <duration index="1" fixed_duration="true">86400</duration>
        //                  </phase_durations>
        //                  <parameters>
        //                    <unlysed_fluid_change_rate units="1/min">1.11667e-2</unlysed_fluid_change_rate>
        //                    <lysed_fluid_change_rate units="1/min">8.33333e-4</lysed_fluid_change_rate>
        //                    <cytoplasmic_biomass_change_rate units="1/min">5.33333e-5</cytoplasmic_biomass_change_rate>
        //                    <nuclear_biomass_change_rate units="1/min">2.16667e-3</nuclear_biomass_change_rate>
        //                    <calcification_rate units="1/min">0</calcification_rate>
        //                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
        //                  </parameters>
        //                </model>
        //              </death>
    }

    private void readVolume(Element el, Phenotype p)
    {
        Volume volume = p.volume;
        for( Element paramElement : getAllElements( el ) )
        {
            switch( paramElement.getTagName() )
            {
                case "total":
                    volume.total = getDoubleVal( paramElement );
                    //                    String total_units = getAttr( paramElement, "units" );
                    break;
                case "fluid_fraction":
                    volume.fluid_fraction = getDoubleVal( paramElement );
                    //                    String fluid_units = getAttr( paramElement, "units" );
                    break;
                case "nuclear":
                    volume.nuclear = getDoubleVal( paramElement );
                    //                    String nuclear_units = getAttr( paramElement, "units" );
                    break;
                case "fluid_change_rate":
                    volume.fluid_change_rate = getDoubleVal( paramElement );
                    //                    String fluid_change_rate_units = getAttr( paramElement, "units" );
                    break;
                case "cytoplasmic_biomass_change_rate":
                    volume.cytoplasmic_biomass_change_rate = getDoubleVal( paramElement );
                    //                    String cytoplasmic_biomass_change_rate_units = getAttr( paramElement, "units" );
                    break;
                case "nuclear_biomass_change_rate":
                    volume.nuclear_biomass_change_rate = getDoubleVal( paramElement );
                    //                    String nuclear_biomass_change_rate_units = getAttr( paramElement, "units" );
                    break;
                case "calcified_fraction":
                    volume.calcified_fraction = getDoubleVal( paramElement );
                    //                    String calcified_fraction_units = getAttr( paramElement, "units" );
                    break;
                case "calcification_rate":
                    volume.calcification_rate = getDoubleVal( paramElement );
                    //                    String calcification_rate_units = getAttr( paramElement, "units" );
                    break;
                case "relative_rupture_volume":
                    volume.relative_rupture_volume = getDoubleVal( paramElement );
                    //                    String relative_rupture_volume_units = getAttr( paramElement, "units" );
                    break;
            }
        }
        volume.adjust();
        p.geometry.update( null, p, 0.0 );
        //              <volume>
        //                <total units="micron^3">2494</total>
        //                <fluid_fraction units="dimensionless">0.75</fluid_fraction>
        //                <nuclear units="micron^3">540</nuclear>
        //                <fluid_change_rate units="1/min">0.05</fluid_change_rate>
        //                <cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
        //                <nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
        //                <calcified_fraction units="dimensionless">0</calcified_fraction>
        //                <calcification_rate units="1/min">0</calcification_rate>
        //                <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
        //              </volume>
    }

    private void readMechanics(Element el, Phenotype p) throws Exception
    {
        Mechanics mechanics = p.mechanics;
        for( Element paramElement : getAllElements( el ) )
        {
            switch( paramElement.getTagName() )
            {
                case "cell_cell_adhesion_strength":
                    mechanics.cell_cell_adhesion_strength = getDoubleVal( paramElement );
                    //                    String cell_cell_adhesion_strength_units = getAttr( paramElement, "units" );
                    break;
                case "cell_cell_repulsion_strength":
                    mechanics.cell_cell_repulsion_strength = getDoubleVal( paramElement );
                    //                    String cell_cell_repulsion_strength_units = getAttr( paramElement, "units" );
                    break;
                case "relative_maximum_adhesion_distance":
                    mechanics.relative_maximum_adhesion_distance = getDoubleVal( paramElement );
                    //                    String relative_maximum_adhesion_distance_units = getAttr( paramElement, "units" );
                    break;
                case "cell_BM_adhesion_strength":
                    mechanics.cell_BM_adhesion_strength = getDoubleVal( paramElement );
                    //                    String cell_BM_adhesion_strength_units = getAttr( paramElement, "units" );
                    break;
                case "cell_BM_repulsion_strength":
                    mechanics.cell_BM_repulsion_strength = getDoubleVal( paramElement );
                    //                    String cell_BM_repulsion_strength_units = getAttr( paramElement, "units" );
                    break;
                case "attachment_elastic_constant":
                    mechanics.attachment_elastic_constant = getDoubleVal( paramElement );
                    //                    String attachment_elastic_constant_units = getAttr( paramElement, "units" );
                    break;
                case "attachment_rate":
                    mechanics.attachment_rate = getDoubleVal( paramElement );
                    //                    String attachment_rate_units = getAttr( paramElement, "units" );
                    break;
                case "detachment_rate":
                    mechanics.detachment_rate = getDoubleVal( paramElement );
                    //                    String detachment_rate_units = getAttr( paramElement, "units" );
                    break;
                case "cell_adhesion_affinities":
                    for( Element adhesionElement : findAllElements( paramElement, "cell_adhesion_affinity" ) )
                    {
                        String target = getAttr( adhesionElement, "name" );
                        double value = getDoubleVal( adhesionElement );
                        int ind = CellDefinition.findCellDefinitionIndex( target );
                        if( ind > -1 )
                            mechanics.cell_adhesion_affinities[ind] = value;
                        else
                            throw new Exception( "Unknown Cell Definition " + target );
                        //                                                { std::cout << "what?!?" << std::endl; }
                    }
                    break;
                case "options":
                    for( Element optionElement : getAllElements( paramElement ) )
                    {
                        switch( optionElement.getTagName() )
                        {
                            case "set_relative_equilibrium_distance":
                                //                                String set_relative_equilibrium_distance_units = getAttr( optionElement, "units" );
                                boolean set_relative_equilibrium_distance_units_isEnabled = getBoolAttr( optionElement, "enabled" );
                                double set_relative_equilibrium_distance_units_value = getDoubleVal( optionElement );
                                if( set_relative_equilibrium_distance_units_isEnabled )
                                    mechanics.set_relative_equilibrium_distance( set_relative_equilibrium_distance_units_value );
                            case "set_absolute_equilibrium_distance":
                                //                                String set_absolute_equilibrium_distance_units = getAttr( optionElement, "units" );
                                boolean set_absolute_equilibrium_distance_units_isEnabled = getBoolAttr( optionElement, "enabled" );
                                double set_absolute_equilibrium_distance_value = getDoubleVal( optionElement );
                                if( set_absolute_equilibrium_distance_units_isEnabled )
                                    mechanics.set_absolute_equilibrium_distance( p, set_absolute_equilibrium_distance_value );
                        }

                    }
            }
        }
        //              <mechanics>
        //                <cell_cell_adhesion_strength units="micron/min">0.4</cell_cell_adhesion_strength>
        //                <cell_cell_repulsion_strength units="micron/min">10.0</cell_cell_repulsion_strength>
        //                <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
        //                <cell_adhesion_affinities>
        //                    <cell_adhesion_affinity name="director cell">1.0</cell_adhesion_affinity>
        //                    <cell_adhesion_affinity name="cargo cell">1.0</cell_adhesion_affinity>
        //                    <cell_adhesion_affinity name="worker cell">1.0</cell_adhesion_affinity>
        //                    </cell_adhesion_affinities>
        //                <options>
        //                  <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
        //                  <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
        //                </options>
        //                <cell_BM_adhesion_strength units="micron/min">4.0</cell_BM_adhesion_strength>
        //                <cell_BM_repulsion_strength units="micron/min">10.0</cell_BM_repulsion_strength>
        //                <attachment_elastic_constant units="1/min">0.5</attachment_elastic_constant>
        //                <attachment_rate units="1/min">10.0</attachment_rate>
        //                <detachment_rate units="1/min">0.0</detachment_rate>
        //              </mechanics
    }

    private void readMotility(Element el, CellDefinition cd, Microenvironment m)
    {
        Phenotype p = cd.phenotype;
        Motility motility = p.motility;

        for( Element child : getAllElements( el ) )
        {
            switch( child.getTagName() )
            {
                case "speed":
                    motility.migration_speed = getDoubleVal( child );
                    //                    String speed_units = getAttr( child, "units" );
                    break;
                case "persistence_time":
                    motility.persistence_time = getDoubleVal( child );
                    //                    String persistence_time_units = getAttr( child, "units" );
                    break;
                case "migration_bias":
                    motility.migration_bias = getDoubleVal( child );
                    //                    String migration_bias_units = getAttr( child, "units" );
                    break;
                case "options":
                {
                    motility.is_motile = getBoolVal( findElement( child, "enabled" ) );
                    motility.restrict_to_2D = getBoolVal( findElement( child, "use_2D" ) );
                    if( m.options.simulate_2D && !motility.restrict_to_2D )
                    {
                        System.out.println( "Note: Overriding to set cell motility for " + cd.name
                                + " to 2D based on microenvironment domain settings ... " );
                        motility.restrict_to_2D = true;
                        break;
                    }
                    Element chemotaxisElement = findElement( child, "chemotaxis" );
                    if( chemotaxisElement != null )
                    {
                        if( getBoolVal( findElement( chemotaxisElement, "enabled" ) ) )// enabled? if so, set the standard chemotaxis function
                            cd.functions.update_migration_bias = new StandardModels.chemotaxis_function();

                        // search for the right chemo index               
                        String substrate_name = getVal( findElement( chemotaxisElement, "substrate" ) );

                        motility.chemotaxis_index = m.findDensityIndex( substrate_name );
                        if( motility.chemotaxis_index < 0 )
                        {
                            System.out.println( "Error: parsing phenotype:motility:options:chemotaxis:  invalid substrate" );
                            System.out.println( "Substrate " + substrate_name + " was not found in the microenvironment." );
                        }
                        String actual_name = m.density_names[motility.chemotaxis_index];
                        if( !substrate_name.equals( actual_name ) )
                        {
                            System.out.println( "Error: attempted to set chemotaxis to \"" + substrate_name
                                    + "\", which was not found in the microenvironment."
                                    + " Please double-check your substrate name in the config file." );
                        }
                        motility.chemotaxis_direction = getIntVal( findElement( chemotaxisElement, "direction" ) );
                    }
                    Element advancedChemotaxisElement = findElement( child, "advanced_chemotaxis" );
                    if( advancedChemotaxisElement != null )
                    {
                        if( getBoolVal( findElement( advancedChemotaxisElement, "enabled" ) ) )
                        {
                            //                            cd.functions.update_migration_bias = new StandardModels.chemotaxis_function();
                            if( cd.functions.update_migration_bias instanceof StandardModels.chemotaxis_function )
                            {
                                System.out.println( "Warning: when processing motility for " + cd.name + " cells: \n"
                                        + "\tBoth chemotaxis and advanced_chemotaxis are enabled.\n"
                                        + "\tThe settings for advanced_chemotaxis override those of chemotaxis." );
                            }

                            if( getBoolAttr( advancedChemotaxisElement, "normalize_each_gradient" ) )
                                cd.functions.update_migration_bias = new StandardModels.advanced_chemotaxis_function_normalized();
                            else
                                cd.functions.update_migration_bias = new StandardModels.advanced_chemotaxis_function();
                        }
                        Element sensitivityEl = findElement( advancedChemotaxisElement, "chemotactic_sensitivities" );
                        if( sensitivityEl != null )
                        {
                            for( Element sensEl : findAllElements( sensitivityEl, "chemotactic_sensitivity" ) )
                            {
                                String substrate_name = getAttr( sensEl, "substrate" );
                                int index = m.findDensityIndex( substrate_name );
                                String actual_name = "";
                                if( index > -1 )
                                {
                                    actual_name = m.density_names[index];
                                }

                                // error check 
                                if( !substrate_name.equals( actual_name ) )
                                {
                                    System.out.println( "Warning: when processing advanced chemotaxis for " + cd.name + " cells: "
                                            + "\tInvalid substrate " + substrate_name + " specified."
                                            + "\tIgnoring this invalid substrate in the chemotaxis function .. " );
                                }
                                else
                                {
                                    cd.phenotype.motility.chemotactic_sensitivities[index] = getDoubleVal( sensEl );
                                }
                            }
                        }
                        else
                        {
                            System.out.println( "Warning: when processing motility for " + cd.name + " cells: "
                                    + "\tAdvanced chemotaxis requries chemotactic_sensitivities."
                                    + "\tBut you have none. Your migration bias will be the zero vector." );
                        }
                    }
                }
            }
        }
        // display summary for diagnostic help 
        if( cd.functions.update_migration_bias instanceof StandardModels.chemotaxis_function && motility.is_motile )
        {
            System.out.println( "Cells of type " + cd.name + " use standard chemotaxis: \n" + "\t d_bias (before normalization) = "
                    + motility.chemotaxis_direction + " * grad(" + m.density_names[motility.chemotaxis_index] + ")" );
        }

        if( cd.functions.update_migration_bias instanceof StandardModels.advanced_chemotaxis_function && motility.is_motile )
        {
            int number_of_substrates = m.density_names.length;

            System.out.println( "Cells of type " + cd.name + " use advanced chemotaxis: \n" + "\t d_bias (before normalization) = "
                    + motility.chemotactic_sensitivities[0] + " * grad(" + m.density_names[0] + ")" );

            for( int n = 1; n < number_of_substrates; n++ )
                System.out.println( motility.chemotactic_sensitivities[n] + " * grad(" + m.density_names[n] + ")" );
        }

        if( cd.functions.update_migration_bias instanceof StandardModels.advanced_chemotaxis_function_normalized && motility.is_motile )
        {
            int number_of_substrates = m.density_names.length;

            System.out.println( "Cells of type " + cd.name + " use normalized advanced chemotaxis: \n"
                    + "\t d_bias (before normalization) = " + motility.chemotactic_sensitivities[0] + " * grad(" + m.density_names[0] + ")"
                    + " / ||grad(" + m.density_names[0] + ")||" );

            for( int n = 1; n < number_of_substrates; n++ )
            {
                System.out.println( motility.chemotactic_sensitivities[n] + " * grad(" + m.density_names[n] + ")" + " / ||grad("
                        + m.density_names[n] + ")||" );
            }
        }
        //              <motility>
        //                <speed units="micron/min">1</speed>
        //                <persistence_time units="min">1</persistence_time>
        //                <migration_bias units="dimensionless">.5</migration_bias>
        //                <options>
        //                  <enabled>false</enabled>
        //                  <use_2D>true</use_2D>
        //                  <chemotaxis>
        //                    <enabled>false</enabled>
        //                    <substrate>director signal</substrate>
        //                    <direction>1</direction>
        //                  </chemotaxis>
        //                  <advanced_chemotaxis>
        //                    <enabled>false</enabled>
        //                    <normalize_each_gradient>false</normalize_each_gradient>
        //                    <chemotactic_sensitivities>
        //                      <chemotactic_sensitivity substrate="director signal">0.0</chemotactic_sensitivity>
        //                      <chemotactic_sensitivity substrate="cargo signal">0.0</chemotactic_sensitivity>
        //                      </chemotactic_sensitivities>
        //                    </advanced_chemotaxis>
        //                </options>
        //              </motility>
    }

    private void readSecretion(Element el, Phenotype p, Microenvironment m)
    {
        Secretion secretion = p.secretion;
        for( Element substrateElement : findAllElements( el, "substrate" ) )
        {
            String name = getAttr( substrateElement, "name" );
            int index = m.findDensityIndex( name );
            String actual_name = m.density_names[index];
            if( !name.equals( actual_name ) )
            {
                System.out.println( "Error: attempted to set secretion/uptake/export for \"" + name
                        + "\", which was not found in the microenvironment.\n"
                        + "       Please double-check your substrate name in the config file." );
            }
            for( Element child : getAllElements( substrateElement ) )
            {
                switch( child.getTagName() )
                {
                    case "secretion_rate":
                        secretion.secretionRates[index] = getDoubleVal( child );
                        //                        String secreion_rate_units = getAttr( child, "units" );
                        break;
                    case "secretion_target":
                        secretion.saturationDensities[index] = getDoubleVal( child );
                        //                        String secretion_target_units = getAttr( child, "units" );
                        break;
                    case "uptake_rate":
                        secretion.uptakeRates[index] = getDoubleVal( child );
                        //                        String uptake_rate_units = getAttr( child, "units" );
                        break;
                    case "net_export_rate":
                        secretion.netExportRates[index] = getDoubleVal( child );
                        //                        String net_export_rate_units = getAttr( child, "units" );
                        break;
                }
            }
        }
        //            // net export rate 
        //            node_sec1 = node_sec.child( "net_export_rate" ); 
        //            if( node_sec1 )
        //            { pS->net_export_rates[index] = xml_get_my_double_value( node_sec1 ); }
        //              <secretion>
        //                <substrate name="director signal">
        //                  <secretion_rate units="1/min">9.9</secretion_rate>
        //                  <secretion_target units="substrate density">1</secretion_target>
        //                  <uptake_rate units="1/min">0</uptake_rate>
        //                  <net_export_rate units="total substrate/min">0</net_export_rate>
        //                </substrate>
        //                <substrate name="cargo signal">
        //                  <secretion_rate units="1/min">0.0</secretion_rate>
        //                  <secretion_target units="substrate density">1</secretion_target>
        //                  <uptake_rate units="1/min">0.0</uptake_rate>
        //                  <net_export_rate units="total substrate/min">0.0</net_export_rate>
        //                </substrate>
        //                </secretion>
    }

    private void readCellInteractions(Element el, CellDefinition cd)
    {
        CellInteractions cellInteractions = cd.phenotype.cell_interactions;
        Element dprElement = findElement( el, "dead_phagocytosis_rate" );
        if( dprElement != null )
        {
            cellInteractions.deadPhagocytosisRate = getDoubleVal( dprElement );
        }
        Element lprsElement = findElement( el, "live_phagocytosis_rates" );
        if( lprsElement != null )
        {
            for( Element lprElement : findAllElements( lprsElement, "phagocytosis_rate" ) )
            {
                // get the name of the target cell type
                String target_name = getAttr( lprElement, "name" );
                // now find its index 
                int index = CellDefinition.findCellDefinitionIndex( target_name );
                // safety first! 
                if( index >= 0 )
                {
                    // if the target is found, set the appropriate rate 
                    cellInteractions.livePhagocytosisRates[index] = getDoubleVal( lprElement );
                    //                    String units = getAttr( lprElement, "units" );
                }
                else
                {
                    System.out.println( "Warning: When processing the " + cd.name + " cell definition: \n" + "\tCould not find cell type "
                            + target_name + " for phagocytosis.\n" + "\tIgnoring this live phagocytosis rate!" );
                }
            }
        }
        Element attackRatesElement = findElement( el, "attack_rates" );
        if( attackRatesElement != null )
        {
            for( Element arElement : findAllElements( attackRatesElement, "attack_rate" ) )
            {
                // get the name of the target cell type
                String target_name = getAttr( arElement, "name" );
                //            // now find its index 
                int index = CellDefinition.findCellDefinitionIndex( target_name );
                if( index >= 0 )
                {
                    cellInteractions.attackRates[index] = getDoubleVal( arElement );
                }
                else
                {
                    System.out.println( "Warning: When processing the " + cd.name + " cell definition: \n" + "\tCould not find cell type "
                            + target_name + " for cell attack.\n" + "\tIgnoring this cell attack rate!" );
                }
            }
        }
        Element fusionRatesElement = findElement( el, "fusion_rates" );
        if( fusionRatesElement != null )
        {
            for( Element drElement : findAllElements( fusionRatesElement, "fusion_rate" ) )
            {
                // get the name of the target cell type
                String target_name = getAttr( drElement, "name" );
                //            // now find its index 
                int index = CellDefinition.findCellDefinitionIndex( target_name );
                if( index >= 0 )
                {
                    cellInteractions.fusionRates[index] = getDoubleVal( drElement );
                }
                else
                {
                    System.out.println( "Warning: When processing the " + cd.name + " cell definition: \n" + "\tCould not find cell type "
                            + target_name + " for cell fusion.\n" + "\tIgnoring this cell fusion rate!" );
                }
            }

            Element damageRateElement = findElement( el, "damage_rate" );
            if( damageRateElement != null )
            {
                cellInteractions.damageRate = getDoubleVal( damageRateElement );
                //                String units = getAttr( damageRateElement, "units" );
            }
        }
        //              <cell_interactions>
        //                <dead_phagocytosis_rate units="1/min">0</dead_phagocytosis_rate>
        //                <live_phagocytosis_rates>
        //                    <phagocytosis_rate name="director cell" units="1/min">0.0</phagocytosis_rate>
        //                    <phagocytosis_rate name="cargo cell" units="1/min">0.0</phagocytosis_rate>
        //                    <phagocytosis_rate name="worker cell" units="1/min">0.0</phagocytosis_rate>
        //                    </live_phagocytosis_rates>
        //
        //                <attack_rates>
        //                      <attack_rate name="director cell" units="1/min">0.0</attack_rate>
        //                      <attack_rate name="cargo cell" units="1/min">0.0</attack_rate>
        //                      <attack_rate name="worker cell" units="1/min">0.0</attack_rate>
        //                      </attack_rates>
        //
        //                <damage_rate units="1/min">1</damage_rate>
        //                <fusion_rates>
        //                      <fusion_rate name="director cell" units="1/min">0.0</fusion_rate>
        //                      <fusion_rate name="cargo cell" units="1/min">0.0</fusion_rate>
        //                      <fusion_rate name="worker cell" units="1/min">0.0</fusion_rate>
        //                      </fusion_rates>
        //
        //              </cell_interactions>
    }

    private void readCellTransformations(Element el, CellDefinition cd)
    {
        CellTransformations transformations = cd.phenotype.cell_transformations;
        Element ratesElement = findElement( el, "transformation_rates" );
        for( Element rateElement : findAllElements( ratesElement, "transformation_rate" ) )
        {
            String name = getAttr( rateElement, "name" );
            int index = CellDefinition.findCellDefinitionIndex( name );
            if( index >= 0 )
            {
                double rate = getDoubleVal( rateElement );
                //              String units = getAttr( rateElement, "units" );
                if( name == cd.name && rate > 1e-16 )
                {
                    System.out.println( "Warning: When processing the " + cd.name + " cell definition:\n" + "\tTransformation from "
                            + cd.name + " to " + name + " is not allowed.\n" + "\tIgnoring this cell transformation rate!" );
                }
                else
                    transformations.transformation_rates[index] = rate;
            }
            else
            {
                System.out.println( "Warning: When processing the " + cd.name + " cell definition:\n" + "\tCould not find cell type " + name
                        + " for cell transformation.\n" + "\tIgnoring this cell transformation rate!" );

            }
        }
        //              <cell_transformations>
        //                <transformation_rates>
        //                    <transformation_rate name="director cell" units="1/min">0.0</transformation_rate>
        //                    <transformation_rate name="cargo cell" units="1/min">0.0</transformation_rate>
        //                    <transformation_rate name="worker cell" units="1/min">0.0</transformation_rate>
        //                    </transformation_rates>
        //                </cell_transformations>
    }

    private void readCustomData(Element el, CellDefinition cd) throws Exception
    {
        for( Element child : getAllElements( el ) )
        {
            CellDefinition defaults = StandardModels.getDefaultCellDefinition();
            // name of teh custom data 
            String name = child.getTagName();

            // units 
            String units = getAttr( child, "units" );
            double[] values; // ( length, 0.0 ); 

            // conserved quantity 
            boolean conserved = getBoolAttr( child, "conserved" );

            // get value(s)
            String str_values = getVal( child );
            values = VectorUtil.csv_to_vector( str_values );

            // add variable if cell defaults  
            // if the custom data is not yet found, add it 
            // first, try scalar 
            if( values.length == 1 )
            {
                // find the variable 
                int n = cd.custom_data.find_variable_index( name );
                // if it exists, overwrite 
                if( n > -1 )
                {
                    cd.custom_data.variables.get( n ).value = values[0];
                }
                // otherwise, add 
                else
                {
                    cd.custom_data.add_variable( name, units, values[0] );
                    defaults.custom_data.add_variable( name, units, values[0] );
                }

                n = cd.custom_data.find_variable_index( name );
            }
            // or vector 
            else
            {
                // find the variable 
                int n = cd.custom_data.find_vector_variable_index( name );
                // if it exists, overwrite 
                if( n > -1 )
                {
                    cd.custom_data.vectorVariables.get( n ).value = values;
                }
                // otherwise, add 
                else
                {
                    cd.custom_data.add_vector_variable( name, units, values );
                }

                n = cd.custom_data.find_vector_variable_index( name );
            }
        }
    }

    public void readUserParameters(Element el, Model model)
    {
        Element parametersElement = findElement( el, "user_parameters" );
        if( parametersElement != null )
        {
            for( Element element : getAllElements( parametersElement ) )
            {
                model.addParameter( element.getTagName(), getVal( element ) );
            }
        }
    }

    public void readInitialConditions(Element el, Model m) throws Exception
    {
        Element parametersElement = findElement( el, "initial_conditions" );
        if( parametersElement != null )
        {
            Element positionsElement = findElement( parametersElement, "cell_positions" );
            boolean enabled = getBoolAttr( positionsElement, "enabled" );
            if( !enabled )
                return;
            String folder = getVal( findElement( positionsElement, "folder" ) );
            String fileName = getVal( findElement( positionsElement, "fileName" ) );
            String input_filename = folder + "/" + fileName;
            CellCSVReader.load_cells_csv( input_filename, m.getMicroenvironment() );
        }
    }

    private Element findElement(Element parent, String tag)
    {
        return findElement( parent.getChildNodes(), tag );
    }

    private Element findElement(NodeList list, String tag)
    {
        for( int i = 0; i < list.getLength(); i++ )
        {
            Node node = list.item( i );
            if( node instanceof Element && ( (Element)node ).getTagName().equals( tag ) )
            {
                return (Element)node;
            }
        }
        return null;
    }

    private List<Element> getAllElements(Element parent)
    {
        List<Element> result = new ArrayList<>();
        NodeList list = parent.getChildNodes();
        for( int i = 0; i < list.getLength(); i++ )
        {
            Node node = list.item( i );
            if( node instanceof Element )
            {
                result.add( (Element)node );
            }
        }
        return result;
    }

    private List<Element> findAllElements(Element parent, String tag)
    {
        List<Element> result = new ArrayList<>();
        for( int i = 0; i < parent.getChildNodes().getLength(); i++ )
        {
            Node node = parent.getChildNodes().item( i );
            if( node instanceof Element && ( (Element)node ).getTagName().equals( tag ) )
            {
                result.add( (Element)node );
            }
        }
        return result;
    }


    private boolean getBoolVal(Element el)
    {
        return Boolean.parseBoolean( getVal( el ) );
    }

    private String getVal(Element el)
    {
        return el.getChildNodes().item( 0 ).getNodeValue();// el.getNodeValue();
    }

    private double getDoubleVal(Element el)
    {
        return Double.parseDouble( getVal( el ) );
    }

    private int getIntVal(Element el)
    {
        return Integer.parseInt( getVal( el ) );
    }

    private double getDoubleAttr(Element el, String name)
    {
        return Double.parseDouble( el.getAttribute( name ) );
    }

    private boolean getBoolAttr(Element el, String name)
    {
        return Boolean.parseBoolean( el.getAttribute( name ) );
    }

    private String getAttr(Element el, String name)
    {
        return el.getAttribute( name );
    }

    private Integer getIntAttr(Element el, String name)
    {
        return Integer.parseInt( el.getAttribute( name ) );
    }

    public static Color readColor(String str)
    {
        switch( str )
        {
            case "green":
                return Color.green;
            case "red":
                return Color.red;
            case "blue":
                return Color.blue;
            case "limegreen":
                return new Color( 50, 205, 50 );
            default:
                return Color.white;
        }

    }
}
