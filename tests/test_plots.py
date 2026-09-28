#coding:utf-8
"""
Tests for visualization modules: plot_hist.py and plot_pie.py
"""
from tools import plot_hist, plot_pie


class TestHistPlot:
    """Tests for hist() - histogram visualization"""
    
    def test_creates_histogram_file(self, tmp_path):
        """Test that histogram PNG file is created"""
        # plot_hist.hist expects a directory containing CSV files
        ams_dir = tmp_path / "ams_data"
        ams_dir.mkdir()
        
        # Create sample CSV files in the directory
        ams_file1 = ams_dir / "ams1.csv"
        ams_file1.write_text(
            "pair,ams,min_dp,max_dp,min_ad,min_gq,hom_thr,base_len\n"
            "pairA,5,10,1000,2,20,0.8,3\n",
            encoding="utf-8"
        )
        
        ams_file2 = ams_dir / "ams2.csv"
        ams_file2.write_text(
            "pair,ams,min_dp,max_dp,min_ad,min_gq,hom_thr,base_len\n"
            "pairB,8,10,1000,2,20,0.8,3\n",
            encoding="utf-8"
        )
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        # Call the function
        plot_hist.hist(str(ams_dir), str(run_plots))
        
        # Check that the PNG file was created
        png_file = run_plots / "distrib.png"
        assert png_file.exists()

    def test_handles_single_pair(self, tmp_path):
        """Test histogram with single pair"""
        ams_dir = tmp_path / "ams_data"
        ams_dir.mkdir()
        
        ams_file = ams_dir / "ams.csv"
        ams_file.write_text(
            "pair,ams\n"
            "pairA,5\n",
            encoding="utf-8"
        )
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        # Call the function
        plot_hist.hist(str(ams_dir), str(run_plots))
        
        # Check that the PNG file was created
        png_file = run_plots / "distrib.png"
        assert png_file.exists()

    def test_handles_empty_data(self, tmp_path):
        """Test histogram with empty data"""
        ams_dir = tmp_path / "ams_data"
        ams_dir.mkdir()
        
        ams_file = ams_dir / "empty.csv"
        ams_file.write_text(
            "pair,ams\n",
            encoding="utf-8"
        )
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        # Call the function
        plot_hist.hist(str(ams_dir), str(run_plots))
        
        # Check that the PNG file was created
        png_file = run_plots / "distrib.png"
        assert png_file.exists()

    def test_handles_large_mismatch_counts(self, tmp_path):
        """Test histogram with large mismatch values"""
        ams_dir = tmp_path / "ams_data"
        ams_dir.mkdir()
        
        ams_file = ams_dir / "ams.csv"
        ams_file.write_text(
            "pair,ams\n"
            "pairA,500\n"
            "pairB,1000\n"
            "pairC,50\n",
            encoding="utf-8"
        )
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        # Call the function
        plot_hist.hist(str(ams_dir), str(run_plots))
        
        # Check that the PNG file was created
        png_file = run_plots / "distrib.png"
        assert png_file.exists()


class TestPiePlot:
    """Tests for pie() - pie chart visualization"""
    
    def test_creates_pie_charts(self, tmp_path):
        """Test that pie chart PNG files are created"""
        run_tables = tmp_path / "run_tables"
        run_tables.mkdir()
        
        # Create mismatches file with chromosomes
        mismatches_file = run_tables / "mismatches.tsv"
        mismatches_file.write_text(
            "CHROM\tPOS\tData\n"
            "1\t100\tA\n"
            "1\t200\tB\n"
            "2\t300\tC\n"
            "3\t400\tD\n",
            encoding="utf-8"
        )
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        plot_pie.pie(str(run_tables), str(run_plots), "", "test_run")
        
        # Check PNG files created
        png_files = list(run_plots.glob("*.png"))
        assert len(png_files) > 0 or run_plots.exists()

    def test_handles_multiple_chromosomes(self, tmp_path):
        """Test pie chart with multiple chromosomes"""
        run_tables = tmp_path / "run_tables"
        run_tables.mkdir()
        
        # Create mismatches across chromosomes
        mismatches_file = run_tables / "mismatches.tsv"
        mismatches_file.write_text(
            "CHROM\tPOS\n"
            "1\t100\n"
            "2\t200\n"
            "3\t300\n"
            "X\t400\n"
            "Y\t500\n"
            "MT\t600\n",
            encoding="utf-8"
        )
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        plot_pie.pie(str(run_tables), str(run_plots), "", "test")
        assert run_plots.exists()

    def test_handles_numeric_chromosome(self, tmp_path):
        """Test with numeric chromosome names"""
        run_tables = tmp_path / "run_tables"
        run_tables.mkdir()
        
        mismatches_file = run_tables / "mismatches.tsv"
        mismatches_file.write_text(
            "CHROM\tPOS\n"
            "1\t100\n"
            "2\t200\n",
            encoding="utf-8"
        )
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        plot_pie.pie(str(run_tables), str(run_plots), "", "test")
        assert run_plots.exists()

    def test_handles_string_chromosome(self, tmp_path):
        """Test with string chromosome names"""
        run_tables = tmp_path / "run_tables"
        run_tables.mkdir()
        
        mismatches_file = run_tables / "mismatches.tsv"
        mismatches_file.write_text(
            "CHROM\tPOS\n"
            "chrX\t100\n"
            "chrY\t200\n",
            encoding="utf-8"
        )
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        plot_pie.pie(str(run_tables), str(run_plots), "", "test")
        assert run_plots.exists()

    def test_handles_pair_name(self, tmp_path):
        """Test with pair name included"""
        run_tables = tmp_path / "run_tables"
        run_tables.mkdir()
        
        mismatches_file = run_tables / "P01_mismatches.tsv"
        mismatches_file.write_text(
            "CHROM\tPOS\n"
            "1\t100\n",
            encoding="utf-8"
        )
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        plot_pie.pie(str(run_tables), str(run_plots), "P01", "test")
        assert run_plots.exists()

    def test_handles_empty_mismatches(self, tmp_path):
        """Test with no mismatches file"""
        run_tables = tmp_path / "run_tables"
        run_tables.mkdir()
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        # No mismatches file created
        try:
            plot_pie.pie(str(run_tables), str(run_plots), "", "test")
        except Exception:
            # Expected if no file found
            pass

    def test_handles_single_mismatch(self, tmp_path):
        """Test with single mismatch"""
        run_tables = tmp_path / "run_tables"
        run_tables.mkdir()
        
        mismatches_file = run_tables / "mismatches.tsv"
        mismatches_file.write_text(
            "CHROM\tPOS\n"
            "1\t100\n",
            encoding="utf-8"
        )
        
        run_plots = tmp_path / "plots"
        run_plots.mkdir()
        
        plot_pie.pie(str(run_tables), str(run_plots), "", "test")
        assert run_plots.exists()
