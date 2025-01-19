# Cargar las librerias ------------------------------------------
import os
import datetime
from SCPipeline.utils import setup_logging
from SCPipeline.preprocessingV3 import SCData, SCDataConfig
def main():
    

    # Obtener las rutas de los directorios desde las variables de entorno ------
    log_dir = os.getenv("LOG_DIR", "./logs")  
    data_dir = os.getenv("DATA_DIR", "./data")  
    output_dir = os.getenv("OUTPUT_DIR", "./data")

    # Asegurarse de que los directorios existen
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    # Usar los directorios en el script
    print(f"Directorio de logs: {log_dir}")
    print(f"Directorio de datos: {data_dir}")
    print(f"Directorio de resultados: {output_dir}")

    # Configuración del logger ----------------------------------------------------
    
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file = os.path.join(log_dir, f"{current_time}_preprocessing.log")
    logger = setup_logging("Logger", log_file)
    if logger:
        logger.info("Logger configurado correctamente.")


    # Configuración de las  muestras --------------------------------  

    # Filtrar muestras

    sample_names = ['AH0319-', 'MH0001-', 'MH0025-', 'MH0032-', 'MH0040-', 'MH0042-', 'MH0043-T-', 'MH0114-T3-', 'MH0125-', 'MH0151-', 'MH0163-', 'MH0167-T-', 'PM0360-','B1-MH0131-', 'B1-MH0177-', 'B1-MH4031-','B1-Tum0554-','MH0114-T2-', 'MH0126-', 'MH0135-', 'SH0106-', 'AH0308-', 'MH0031-', 'MH0069-', 'MH0161-', 'MH0176-', 'PM0337-']

    
    all_files = os.listdir(data_dir)
    
    samples = [file for file in all_files if any(name in file for name in sample_names)]



    # Crear objeto SCDataConfig -----------------------------------
    adata_config = SCDataConfig(data_dir=data_dir, samples=samples)
  
    # Crear objeto SCData y generar Adata --------------------------
    adata_obj = SCData(Adata_config=adata_config)

    # Preprocesar -----------------------------------------------
    adata_obj.create_adata_list()
    adata_obj.compute_qc_metrics()
    adata_obj.detect_and_remove_doublets(remove = False) 
    adata_obj.calculate_percentiles()
    adata_obj.calculate_thresholds()
    adata_obj.filter_low_quality_cells()

    modelos = ["Immune_All_High.pkl","Cells_Adult_Breast.pkl"]
    
    adata_obj.annotate_celltypist(modelos_deseados = modelos)

    adata_obj.filter_low_quality_genes()
    adata_obj.detect_and_remove_doublets(remove = True)
    adata_obj.concatenate_samples()
    
    adata_final = adata_obj.adata
    adata_final.write_h5ad(os.path.join(output_dir, "adata.h5ad"))


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error inesperado: {e}")