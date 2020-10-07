#include "plotFile.h"
#include <vector>
#include <arrow/array.h>
#include <arrow/api.h>
#include <memory>
#include <arrow/ipc/writer.h>
#include <parquet/arrow/writer.h>
#include <arrow/io/file.h>
#include "tools.h"


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

void writeSingleLevel(MultiFab & phi_real, MultiFab & phi_imag , Geometry & geom)
{
    std::shared_ptr<arrow::Field> field_phi_real;
    field_phi_real = arrow::field("phi_real", arrow::float64() ) ;
    std::shared_ptr<arrow::Schema> schema;
    schema = arrow::schema({field_phi_real });

    int iB=0;

    arrow::DoubleBuilder field_left_x_builder ;
    arrow::DoubleBuilder field_right_x_builder ;
    arrow::Int32Builder field_index_left_x_builder ;
    arrow::Int32Builder field_index_right_x_builder ;
    arrow::Int32Builder field_box_builder ;
    arrow::Int32Builder field_ghosts_x_builder ;



	const Real* dx = geom.CellSize(); 
	const Real* prob_lo = geom.ProbLo();
    /* write data to boxes */
	for ( MFIter mfi(phi_real); mfi.isValid(); ++mfi ) 
	{ 
	    const Box& bx = mfi.growntilebox();
        const Box& bxValid = mfi.tilebox();
        const int* loValid = bxValid.loVect(); 
	    

	    const int* lo = bx.loVect(); 
	    const int *hi= bx.hiVect(); 
	    Array4< Real> const & phi_real_array = phi_real[mfi].array();
        Array4< Real> const & phi_imag_array = phi_imag[mfi].array();


        size_t  size = hi[0] - lo[0] + 1;
        int j=0;
        int k=0;
        int i=lo[0];

        field_box_builder.Append(iB);
        field_left_x_builder.Append(prob_lo[0] + i * lo[0]*dx[0]);
        field_right_x_builder.Append(prob_lo[0] + i * hi[0]*dx[0]);

        field_index_left_x_builder.Append(lo[0]);
        field_index_right_x_builder.Append(hi[0]);
      
        field_ghosts_x_builder.Append(loValid[0] - lo[0] );




        arrow::DoubleBuilder builder;
        builder.AppendValues( &phi_real_array(i,j,k)     ,  size);

        std::shared_ptr<arrow::Array> phi_real_array_arrow;

        arrow::Status st = builder.Finish(&phi_real_array_arrow);
        auto table = arrow::Table::Make(schema, {phi_real_array_arrow});


        std::shared_ptr<arrow::io::FileOutputStream> outfile;
    PARQUET_ASSIGN_OR_THROW(
      outfile,
      arrow::io::FileOutputStream::Open("out/box=" + std::to_string(iB) ));

    PARQUET_THROW_NOT_OK(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1000000)) ;

      iB++;

    }

    // save box data
    std::shared_ptr<arrow::Array> field_box_array;
    std::shared_ptr<arrow::Array> field_left_x_array;
    std::shared_ptr<arrow::Array> field_right_x_array;
    std::shared_ptr<arrow::Array> field_index_left_x_array;
    std::shared_ptr<arrow::Array> field_index_right_x_array;
    std::shared_ptr<arrow::Array> field_ghosts_x_array;

    field_left_x_builder.Finish(&field_left_x_array);
    field_right_x_builder.Finish(&field_right_x_array);
    field_index_left_x_builder.Finish(&field_index_left_x_array);
    field_index_right_x_builder.Finish(&field_index_right_x_array);
    field_ghosts_x_builder.Finish(&field_ghosts_x_array);
    field_box_builder.Finish(&field_box_array);
    

    



    

    auto boxScheme = createBoxesScheme() ;
    auto table = arrow::Table::Make(boxScheme, {field_box_array, field_left_x_array,field_right_x_array,field_index_left_x_array,field_index_right_x_array,field_ghosts_x_array});

    std::shared_ptr<arrow::io::FileOutputStream> outfile;
    PARQUET_ASSIGN_OR_THROW(
      outfile,
      arrow::io::FileOutputStream::Open("out/boxes" ));

    PARQUET_THROW_NOT_OK(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1000000)) ;






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