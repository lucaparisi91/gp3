#include "plotFile.h"
#include <vector>
#include <arrow/array.h>
#include <arrow/api.h>
#include <memory>
#include <arrow/ipc/writer.h>
#include <parquet/arrow/writer.h>
#include <arrow/io/file.h>
#include "tools.h"
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

auto createBoxesScheme()
{
    auto  field_box = arrow::field("box", arrow::int32() ) ;
    
    auto  field_left_x = arrow::field("left_x", arrow::float64() ) ;
    auto  field_right_x = arrow::field("right_x", arrow::float64() ) ;
    auto field_index_left_x = arrow::field("index_left_x", arrow::int32() );
    auto field_index_right_x = arrow::field("index_right_x", arrow::int32() ) ;
    auto  field_ghosts_x = arrow::field("ghosts_x", arrow::int32() ) ;
    
    return arrow::schema({field_box,field_left_x, field_right_x, field_index_left_x, field_index_right_x , field_ghosts_x});

}

void writeSingleLevel(MultiFab & phi_real, MultiFab & phi_imag , Geometry & geom, const std::string & dirname, Real time)
{
    json_t jBoxes;

    if ( not fs::exists(dirname) )
    {
        fs::create_directory(dirname);
    }



    std::shared_ptr<arrow::Field> field_phi_real;
    field_phi_real = arrow::field("phi_real", arrow::float64() ) ;
    std::shared_ptr<arrow::Schema> schema;
    schema = arrow::schema({field_phi_real });

    int iB=0;
    



	const Real* dx = geom.CellSize(); 
	const Real* prob_lo = geom.ProbLo();
    /* write data to boxes */
	for ( MFIter mfi(phi_real); mfi.isValid(); ++mfi ) 
	{ 
	    const Box& bx = mfi.growntilebox();
        const Box& bxValid = mfi.tilebox();
        const int* loValid = bxValid.loVect();
        const int* hiValid = bxValid.hiVect();

	    const int* lo = bx.loVect(); 
	    const int *hi= bx.hiVect(); 
	    Array4< Real> const & phi_real_array = phi_real[mfi].array();
        Array4< Real> const & phi_imag_array = phi_imag[mfi].array();

        json_t jBox;
        jBox["shape"]= {  hiValid[0] - loValid[0] + 1  } ;
        jBox["lower_index"] ={loValid[0]};
        jBox["ghosts"]= { {  loValid[0] -lo[0] , hi[0] - hiValid[0] }  } ;
        jBox["domain"]={ { prob_lo[0] + loValid[0]*dx[0], prob_lo[0] + (hiValid[0] + 1 )*dx[0]   }};
        jBox["index"] = 0;
        jBoxes.push_back(jBox);


        // save the fab data in parquet files
        size_t  size = hi[0] - lo[0] + 1;
        int j=0;
        int k=0;
        int i=lo[0];

        arrow::DoubleBuilder builder;
        builder.AppendValues( &phi_real_array(i,j,k)     ,  size);

        std::shared_ptr<arrow::Array> phi_real_array_arrow;

        arrow::Status st = builder.Finish(&phi_real_array_arrow);
        auto table = arrow::Table::Make(schema, {phi_real_array_arrow});


        std::shared_ptr<arrow::io::FileOutputStream> outfile;
    PARQUET_ASSIGN_OR_THROW(
      outfile,
      arrow::io::FileOutputStream::Open(dirname + "/box=" + std::to_string(iB) + ".parquet" ));

    PARQUET_THROW_NOT_OK(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1000000)) ;

      iB++;

    }

    // saves json settings
    json_t jDescriptions;
    jDescriptions["boxes"] = jBoxes;
    if (geom.Coord() == 2) 
    {
        jDescriptions["coordinates"]="spherical";
    }
    else if (geom.Coord()==0 )
    {
        jDescriptions["coordinates"]="cartesian";
    }
    jDescriptions["time"]=time;

    std::ofstream f;
    f.open(dirname + "/description.json");
    f << jDescriptions ;

    f.close();

}



void writeSingleLevelTest()
{
    int N=100;
    std::vector<double> data ;

    data.resize( N);
    for (int i=0;i<N;i++)
    {
        data[i] = i * 1./N ;
    }

    arrow::DoubleBuilder builder;
    builder.AppendValues(data.data(),N);

    std::shared_ptr<arrow::Array> dataArray;
    arrow::Status st = builder.Finish(&dataArray);

    std::shared_ptr<arrow::Field> field_x;
    field_x = arrow::field("A", arrow::float64() ) ;
    std::shared_ptr<arrow::Schema> schema;

    schema = arrow::schema({field_x });

    std::shared_ptr<arrow::Table> table; 
    table = arrow::Table::Make(schema, {dataArray});

     
    std::string filename="out.parquet";

    std::shared_ptr<arrow::io::FileOutputStream> outfile;
  PARQUET_ASSIGN_OR_THROW(
      outfile,
      arrow::io::FileOutputStream::Open(filename));


  PARQUET_THROW_NOT_OK(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 3)) ;
     

}