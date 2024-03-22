package ru.biosoft.physicell.fba;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Intracellular;
import ru.biosoft.physicell.core.Phenotype;

public class IntracellularFBA extends Intracellular
{
    public FBA_model model = new FBA_model();

    public Map<String, String> parameterMapping = new HashMap<>();

    Map<String, Double> parameters = new HashMap<>();;
    Map<String, exchange_data> substrate_exchanges = new HashMap<>();

    //    Map<String, String> nameMapping = new HashMap<>();

    double next_model_run = 0;


    public IntracellularFBA()
    {
        this.intracellular_type = "dfba";
    }

    public static class kinetic_parm
    {
        String name;
        String units;
        double value;
    }

    public static class exchange_data
    {
        String density_name;
        String fba_flux_id;
        int density_index;
        kinetic_parm Km;
        kinetic_parm Vmax;
    };

    IntracellularFBA(IntracellularFBA copy)
    {
        intracellular_type = copy.intracellular_type;
        //        sbml_filename = copy.sbml_filename;
        parameters = copy.parameters;
        // model = copy.model;
        //        model.readSBMLModel( copy.sbml_filename.c_str() );
        //        model.initLpModel();
        model.runFBA();
    }

    public String exchange_flux_density_map(String name)
    {
        return parameterMapping.get( name );
    }

    @Override
    public IntracellularFBA clone()
    {
        IntracellularFBA result = new IntracellularFBA();
        result.model = model.clone();
        result.parameterMapping.putAll( parameterMapping );
        result.parameters.putAll( parameters );
        result.substrate_exchanges.putAll( substrate_exchanges );
        result.parameterMapping.putAll( parameterMapping );
        return result;
    }

    //void initialize_intracellular_from_pugixml(pugi::xml_node& node)
    //{
    //
    //    pugi::xml_node node_sbml = node.child( "sbml_filename" );
    //
    //	if ( node_sbml )
    //	{ 
    //        this.sbml_filename = PhysiCell::xml_get_my_string_value (node_sbml);
    //        
    //    }
    //    else
    //       throw new Exception( "Error: attempted get sbml_filename but not foun. Please double-check your exchange nodes in the config file." + std::endl;
    //	
    //    pugi::xml_node node_exchange = node.child( "exchange" );
    //    
    //	while( node_exchange )
    //	{
    //        exchange_data ex_struct;
    //        kinetic_parm Km;
    //        kinetic_parm Vmax;
    //        
    //
    //		string density_name = node_exchange.attribute( "substrate" ).value(); 
    //        int density_index = microenvironment.find_density_index( density_name ); 
    ////        std::cout + "Parsing " + density_name + std::endl;
    //        String actual_name = microenvironment.density_names[ density_index ]; 
    //			
    //        // error check 
    //        if( std::strcmp( density_name.c_str() , actual_name.c_str() ) != 0 )
    //        {
    //            throw new Exception("Error: attempted to set secretion/uptake/export for " 
    //                + density_name + ", which was not found in the microenvironment. Please double-check your substrate name in the config file.");
    //        }
    //        
    //        pugi::xml_node node_fba_flux = node_exchange.child( "fba_flux" ); 
    //		if( node_fba_flux )
    //		{  
    //            ex_struct.fba_flux_id = PhysiCell::xml_get_my_string_value(node_fba_flux);
    //        }
    //        else 
    //            throw new Exception("Error: attempted get fba_flux node for " 
    //             + ex_struct.density_name + ", but not found. Please double-check your exchange nodes in the config file."); 
    //
    //        pugi::xml_node node_Km = node_exchange.child( "Km" ); 
    //		if( node_Km )
    //		{
    //            Km.name = "Km";
    //            Km.untis = node_Km.attribute("units").value();
    //            Km.value = PhysiCell::xml_get_my_double_value(node_Km);
    //            
    //        }
    //        else
    //            throw new Exception("Error: attempted get Km node for "
    //            + density_name + ", but not found. Please double-check your exchange nodes in the config file."); 
    //
    //        pugi::xml_node node_Vmax = node_exchange.child( "Vmax" ); 
    //		if( node_Vmax )
    //		{
    //            Vmax.name = "Vmax";
    //            Vmax.untis = node_Vmax.attribute("units").value();
    //            Vmax.value = PhysiCell::xml_get_my_double_value(node_Vmax);
    //            
    //        }
    //        else {
    //            throw new Exception("Error: attempted get Vmax node for "
    //           + density_name + ", but not found. Please double-check your exchange nodes in the config file.");
    //        }
    //		
    //        ex_struct.density_name = density_name;
    //        ex_struct.density_index = density_index;
    //        ex_struct.Km = Km;
    //        ex_struct.Vmax = Vmax;
    //
    //        this.substrate_exchanges[density_name] = ex_struct;
    //		node_exchange = node_exchange.next_sibling( "exchange" ); 
    //	}
    //
    //   System.out.println( "Loaing SBML model from: " + this.sbml_filename );
    //    this.model.readSBMLModel(this.sbml_filename.c_str());
    //    this.model.initLpModel();
    //    this.model.runFBA();
    //
    //}

    @Override
    public void start()
    {
        // return 0;
    }

    @Override
    public boolean need_update(double curTime)
    {
        return curTime >= this.next_model_run;
    }

    @Override
    public void update(Cell pCell, Phenotype phenotype, double dt)
    {
        for( Entry<String, exchange_data> entry : substrate_exchanges.entrySet() )
        {
            String substrate_name = entry.getKey();
            exchange_data ex_strut = entry.getValue();

            double Vmax = ex_strut.Vmax.value;
            double Km = ex_strut.Km.value;

            // geting the amount of substrate
            double substrate_conc = pCell.nearest_density_vector()[ex_strut.density_index];
            // useing irreversible Michaelis Menten kinetics to estimate the flux bound
            double flux_bound = ( Vmax * substrate_conc ) / ( Km + substrate_conc ); // should be calculated from density
            // Change sign to use as lower bound of the exchange flux
            flux_bound *= -1;
            // Updateing the lower bound of the corresponding exchange flux

            System.out.println( " - [" + substrate_name + "] = " + substrate_conc );
            System.out.println( " ==> " + ex_strut.fba_flux_id + " = " + flux_bound );

            this.model.setReactionLowerBound( ex_strut.fba_flux_id, flux_bound );
        }
        //        this.model.runFBA();
        // this.update_phenotype_parameters(phenotype);//TODO???
        // return 0;
    }

public int update_phenotype_parameters(Phenotype phenotype)
{
    /*
        Required steps for the update

        1- rescale exchange fluxes from the dfba model
        2- update the net_export_rates using the rescaled exchanges
        3- remove the internalized substrates if needed
        4- update the cell volumne using the growth rate from FBA
        5- if no-growth check non-growth associated mantainanace to update apoptotic rate
    */
    for( Entry<String, exchange_data> entry : substrate_exchanges.entrySet() )
    {
        String substrate_name = entry.getKey();
        exchange_data ex_strut = entry.getValue();
        String fba_flux_id = ex_strut.fba_flux_id;
        //        FBA_reaction exchange = this.model.getReaction(fba_flux_id);
        //        double flux = exchange.getFluxValue();
        double flux = model.getFlux( fba_flux_id );
        // how to rescale FBA exchanges into net_export_rates
        float scaling = 1;
        flux *= scaling;
        phenotype.secretion.netExportRates[ex_strut.density_index] = flux;
    }

    double delta_vol = 1;
    double growth_rate = 1; //
    //    phenotype.volume.total; // FBA predicts growth rate mu
    phenotype.volume.multiply_by_ratio(delta_vol);

    return 0;
}



    void save_dFBA(String path, String index)
    {

    }


    @Override
    public void addPhenotypeSpecies(String code, String species)
    {
        // TODO Auto-generated method stub

    }


    @Override
    public void step() throws Exception
    {
        // TODO Auto-generated method stub

    }


    @Override
    public void inherit(Cell cell)
    {
        // TODO Auto-generated method stub

    }


    @Override
    public double getParameterValue(String name) throws Exception
    {
        // TODO Auto-generated method stub
        return 0;
    }


    @Override
    public void setParameterValue(String name, double value) throws Exception
    {
        // TODO Auto-generated method stub

    }


    @Override
    public void setDT(double dt)
    {
        // TODO Auto-generated method stub

    }


    @Override
    public int updatePhenotypeParameters(Microenvironment microenvirionment, Phenotype phenotype) throws Exception
    {
        // TODO Auto-generated method stub
        return 0;
    }


    @Override
    public int validate_PhysiCell_tokens(Microenvironment microenvirionment, Phenotype phenotype) throws Exception
    {
        // TODO Auto-generated method stub
        return 0;
    }

}