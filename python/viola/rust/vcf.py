import pandas as pd
class Vcf:
    def __init__(self, rust_vcf):
        self.rust_vcf = rust_vcf
    
    @property
    def sv_count(self):
        return self.rust_vcf.sv_count()

    @property
    def contigs(self):
        return self.rust_vcf.contigs()
    
    @property
    def ids(self):
        return self.rust_vcf.ids()
    
    def _repr_html_(self):
        return self.view()
    
    def view(self):
        doc_link = 'https://dermasugita.github.io/ViolaDocs/docs/html/reference/vcf.html'
        header = f"""<h1>Vcf Object</h1>
        <h2>Summary</h2>
        <p>Number of SVs: {self.sv_count}</p>
        <a href={doc_link}>Documentation link</a>
        """

        positions = self.get_positions_table()
        head = positions.head()

        head_ids = head["id"]
        head_be1 = head["chrom1"] + ':' + head["pos1"].astype(str)
        head_be2 = head["chrom2"] + ':' + head["pos2"].astype(str)
        head_strand = head["strand1"] + head["strand2"]
        head_qual = head["qual"]
        head_svtype = head["svtype"]
        
        df_ret = pd.DataFrame({
            "id": head_ids,
            "be1": head_be1,
            "be2": head_be2,
            "strand": head_strand,
            "qual": head_qual,
            "svtype": head_svtype
        })
        df_repr = df_ret._repr_html_()
        return header + df_repr
    
    def get_table(self, table_name):
        if table_name == "positions":
            return self.get_positions_table()
    
    def get_positions_table(self):
        positions = self.rust_vcf.get_positions_table()
        columns = positions.get_columns()
        data = {col: getattr(positions, col) for col in columns}
        return pd.DataFrame(data)
        
    def get_str_info_table(self):
        str_info = self.rust_vcf.get_info_tables().flag_info
        columns = str_info.get_columns()
        data = {col: getattr(str_info, col) for col in columns}
        return pd.DataFrame(data)
    
    def get_str_format_table(self):
        str_format = self.rust_vcf.get_format_tables().int_format
        columns = str_format.get_columns()
        data = {col: getattr(str_format, col) for col in columns}
        return pd.DataFrame(data)
        
    
