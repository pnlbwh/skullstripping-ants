#ifndef utilitiesIO_hxx_
#define utilitiesIO_hxx_

#include <csignal>

// itk
#include "itkAffineTransform.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"

#include "itkVector.h"

#include "itkTransformFactoryBase.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"



// vnl
#include "vnl/vnl_matrix.h"


// local
#include "utilitiesIO.h"

namespace kalmanAtlas
{
  /**
   * readImage
   */
  template< typename itkImage_t >
  typename itkImage_t::Pointer readImage(const char *fileName)
  {
    typedef itk::ImageFileReader< itkImage_t > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename itkImage_t::Pointer image;
    
    try
      {
        reader->Update();
        image = reader->GetOutput();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cerr<< "ExceptionObject caught !" << std::endl; 
        std::cerr<< err << std::endl; 
        raise(SIGABRT);
      }

    return image;
  }


  /**
   * writeImage
   */
  template< typename itkImage_t > void writeImage(typename itkImage_t::Pointer img, const char *fileName)
  {
    typedef itk::ImageFileWriter< itkImage_t > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(img);
    writer->UseCompressionOn();

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl; 
        std::cout << err << std::endl; 
        raise(SIGABRT);   
      }
  }


  /**
   * Read a series of images.
   */
  template< typename itkImage_t > 
  std::vector< typename itkImage_t::Pointer >
  readImageSeries( const std::vector< std::string >& imageNameList )
  {
    typedef typename itkImage_t::Pointer itkImagePointer_t;
    typedef std::vector< itkImagePointer_t > itkImageList_t;
    typedef itk::ImageFileReader< itkImage_t > itkImageReader_t;


    itkImageList_t imageSeries;
    
    int n = imageNameList.size();
    for (int i = 0; i < n; ++i)
      {
        std::string thisName = imageNameList[i];

        typename itkImageReader_t::Pointer reader = itkImageReader_t::New();
        reader->SetFileName(thisName);

        itkImagePointer_t img;

        try
          {
            reader->Update();
            img = reader->GetOutput();
          }
        catch ( itk::ExceptionObject &err)
          {
            std::cerr<< "ExceptionObject caught !" << std::endl; 
            std::cerr<< err << std::endl; 
            raise(SIGABRT);
          }


        imageSeries.push_back(img);
      }

    return imageSeries;
  }

  /*
   * readTextLineToListOfString   
   */
  template<typename TNull>
  std::vector< std::string > readTextLineToListOfString(const char* textFileName)
  {
    /* The file MUST end with an empty line, then each line will be
       stored as an element in the returned vector object. */


    // Here speed and space is not a big issue, so just let it copy and return stuff.
    std::vector< std::string > listOfStrings;

    std::ifstream f(textFileName);
    std::string thisLine;

    if (f.good())
      {
        while( std::getline(f, thisLine) )
          {
            listOfStrings.push_back(thisLine);
          }
      }
    else
      {
        std::cerr<<"Error: can not open file:"<<textFileName<<std::endl;
        raise(SIGABRT);
      }

    f.close();

    return listOfStrings;
  }


  /**
   * write a component of a vector image
   */
  template< typename itkVectorImage_t > 
  void 
  writeVectorImage(typename itkVectorImage_t::Pointer img, const char *fileName, int component)
  {
    typedef itk::Image<double, itkVectorImage_t::ImageDimension> itkImage_t;
    typename itkImage_t::Pointer componentImg = itkImage_t::New();
    componentImg->SetRegions(img->GetLargestPossibleRegion() );
    componentImg->Allocate();


    typedef itk::ImageRegionIterator<itkVectorImage_t> VectorIteratorType;
    VectorIteratorType vIter(img, img->GetLargestPossibleRegion());

    typedef itk::ImageRegionIterator<itkImage_t> IteratorType;
    IteratorType iter(componentImg, componentImg->GetLargestPossibleRegion());

    for (vIter.GoToBegin(), iter.GoToBegin(); !vIter.IsAtEnd(); ++iter, ++vIter)
      {
        iter.Set(vIter.Get()[component]);
      }

    typedef itk::ImageFileWriter< itkImage_t > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(componentImg);
    writer->UseCompressionOn();

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl; 
        std::cout << err << std::endl; 
        abort();   
      }
  }


  /**
   * write affine transformation from file
   */
  template< typename CoordRepType, unsigned int NDimensions >
  void writeAffineTransformationToFile(typename itk::AffineTransform<CoordRepType, NDimensions>::Pointer trans, const char *fileName)
  {
    itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
    writer->SetInput(trans);
    writer->SetFileName(fileName);
    writer->Update();

    return;
  }


  /**
   * read affine transformation file
   */
  template< typename CoordRepType, unsigned int NDimensions >
  typename itk::AffineTransform<CoordRepType, NDimensions>::Pointer  
  readAffineTransformFromFile(const char *fileName)
  {
    // Register default transforms
    itk::TransformFactoryBase::RegisterDefaultTransforms();

    typedef itk::TransformFileReader TransformReaderType;
    typedef TransformReaderType::TransformListType TransformListType;
    typedef TransformReaderType::TransformType BaseTransformType;
 
    TransformReaderType::Pointer reader = TransformReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();

    //TransformListType* list = reader->GetTransformList();
    //BaseTransformType* transform = list->front();

    typedef itk::AffineTransform<CoordRepType, NDimensions> TransformType;
    // TransformType* trans = dynamic_cast< TransformType* >( transform );
                                                                                
    typedef itk::TransformFileReader::TransformListType* TransformListPointer;
    TransformListPointer transforms = reader->GetTransformList();                      
    itk::TransformFileReader::TransformListType::const_iterator it = transforms->begin();                                                       
    if (!strcmp((*it)->GetNameOfClass(), "AffineTransform"))                        
      {
        typename TransformType::Pointer affineTransform = static_cast<TransformType*>((*it).GetPointer());

        return affineTransform;
      }

    abort();
  }



}// kalmanAtlas

#endif
