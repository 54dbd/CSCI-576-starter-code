#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Image Error Analysis Script
Analyzes error trends in C1M1, C1M2, C2M1, C2M2 folders

Image naming format:
out_<ColorSpace>_<quantizeMethod>__N<total_bit_budget>_<Q1Q2Q3_bit_for_3_channel>_error_<error_value>.png

Where:
- ColorSpace: RGB, YUV
- quantizeMethod: equal, nonuniform
- total_bit_budget: N4, N6, N8
- Q1Q2Q3: bit allocation for three channels, e.g., 112, 121, 211
- error_value: error value
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set font to avoid character issues
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False

class ImageErrorAnalyzer:
    def __init__(self, base_dir="."):
        self.base_dir = Path(base_dir)
        self.folders = ['C1M1', 'C1M2', 'C2M1', 'C2M2']
        self.data = []
        
    def parse_filename(self, filename):
        """
        Parse image filename to extract parameters
        
        Format: out_<ColorSpace>_<quantizeMethod>__N<total_bit_budget>_<Q1Q2Q3_bit_for_3_channel>_error_<error_value>.png
        """
        # Use regex to parse filename, note the double underscore
        pattern = r'out_(RGB|YUV)_(equal|nonuniform)_N(\d+)_(\d{3})_error_(\d+)\.png'
        match = re.match(pattern, filename)
        
        if match:
            colorspace = match.group(1)
            quantize_method = match.group(2)
            total_bits = int(match.group(3))
            channel_bits = match.group(4)
            error_value = int(match.group(5))
            
            # Parse channel bit allocation
            q1, q2, q3 = int(channel_bits[0]), int(channel_bits[1]), int(channel_bits[2])
            
            return {
                'colorspace': colorspace,
                'quantize_method': quantize_method,
                'total_bits': total_bits,
                'channel_bits': channel_bits,
                'q1': q1,
                'q2': q2,
                'q3': q3,
                'error_value': error_value,
                'filename': filename
            }
        return None
    
    def collect_data(self):
        """Collect data from all folders"""
        print("Collecting data...")
        
        for folder in self.folders:
            folder_path = self.base_dir / folder
            if not folder_path.exists():
                print(f"Warning: Folder {folder} does not exist")
                continue
                
            print(f"Processing folder: {folder}")
            files = list(folder_path.glob("*.png"))
            
            for file in files:
                parsed_data = self.parse_filename(file.name)
                if parsed_data:
                    parsed_data['folder'] = folder
                    self.data.append(parsed_data)
                else:
                    print(f"Cannot parse filename: {file.name}")
        
        print(f"Total collected {len(self.data)} valid data points")
        return self.data
    
    def create_dataframe(self):
        """Create DataFrame for analysis"""
        if not self.data:
            self.collect_data()
        
        if not self.data:
            print("Error: No data collected")
            return pd.DataFrame()
        
        df = pd.DataFrame(self.data)
        
        # Check if required columns exist
        required_columns = ['folder', 'total_bits', 'error_value', 'q1', 'q2', 'q3']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f"Error: Missing required columns: {missing_columns}")
            return pd.DataFrame()
        
        # Add calculated columns
        df['config'] = df['folder']  # Configuration identifier
        df['bit_efficiency'] = df['total_bits'] / df['error_value']  # Bit efficiency
        df['channel_sum'] = df['q1'] + df['q2'] + df['q3']  # Channel bit sum
        
        return df
    
    def analyze_error_trends(self):
        """Analyze error trends"""
        df = self.create_dataframe()
        
        if df.empty:
            print("Error: Cannot create DataFrame, please check data")
            return df
        
        print("\n=== Data Overview ===")
        print(f"Total data points: {len(df)}")
        print(f"Configuration distribution:")
        print(df['config'].value_counts())
        print(f"\nColor space distribution:")
        print(df['colorspace'].value_counts())
        print(f"\nQuantization method distribution:")
        print(df['quantize_method'].value_counts())
        print(f"\nTotal bit budget distribution:")
        print(df['total_bits'].value_counts())
        
        print("\n=== Error Value Statistics ===")
        for config in df['config'].unique():
            config_data = df[df['config'] == config]
            print(f"\n{config}:")
            print(f"  Error range: {config_data['error_value'].min():,} - {config_data['error_value'].max():,}")
            print(f"  Mean Error: {config_data['error_value'].mean():,.0f}")
            print(f"  Median Error: {config_data['error_value'].median():,.0f}")
            print(f"  Std Dev: {config_data['error_value'].std():,.0f}")
        
        return df
    
    def create_visualizations(self, df):
        """Create visualization charts"""
        print("\nGenerating visualization charts...")
        
        if df.empty:
            print("Error: DataFrame is empty, cannot generate charts")
            return
        
        # Set chart style
        plt.style.use('default')
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Image Error Analysis Results', fontsize=16, fontweight='bold')
        
        # 1. Error distribution by configuration
        configs = df['config'].unique()
        error_data = [df[df['config'] == config]['error_value'].values for config in configs]
        axes[0, 0].boxplot(error_data, labels=configs)
        axes[0, 0].set_title('Error Distribution by Configuration')
        axes[0, 0].set_xlabel('Configuration')
        axes[0, 0].set_ylabel('Error Value')
        axes[0, 0].tick_params(axis='x', rotation=45)
        
        # 2. Total bit budget vs Error scatter plot
        for config in df['config'].unique():
            config_data = df[df['config'] == config]
            axes[0, 1].scatter(config_data['total_bits'], config_data['error_value'], 
                              label=config, alpha=0.7, s=50)
        axes[0, 1].set_title('Total Bit Budget vs Error Value')
        axes[0, 1].set_xlabel('Total Bit Budget')
        axes[0, 1].set_ylabel('Error Value')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. Color space comparison
        colorspace_data = df.groupby(['colorspace', 'config'])['error_value'].mean().unstack()
        colorspace_data.plot(kind='bar', ax=axes[1, 0])
        axes[1, 0].set_title('Error Comparison by Color Space')
        axes[1, 0].set_xlabel('Color Space')
        axes[1, 0].set_ylabel('Mean Error Value')
        axes[1, 0].legend(title='Configuration')
        axes[1, 0].tick_params(axis='x', rotation=0)
        
        # 4. Quantization method comparison
        quantize_data = df.groupby(['quantize_method', 'config'])['error_value'].mean().unstack()
        quantize_data.plot(kind='bar', ax=axes[1, 1])
        axes[1, 1].set_title('Error Comparison by Quantization Method')
        axes[1, 1].set_xlabel('Quantization Method')
        axes[1, 1].set_ylabel('Mean Error Value')
        axes[1, 1].legend(title='Configuration')
        axes[1, 1].tick_params(axis='x', rotation=0)
        
        plt.tight_layout()
        plt.savefig('image_error_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("Visualization charts saved as: image_error_analysis.png")
    
    def create_trend_visualizations(self, df):
        """Create trend line charts to show error changes"""
        print("\nGenerating trend line charts...")
        
        if df.empty:
            print("Error: DataFrame is empty, cannot generate trend charts")
            return
        
        # Set chart style
        plt.style.use('default')
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Error Trend Analysis', fontsize=16, fontweight='bold')
        
        # 1. Error trend by total bit budget
        axes[0, 0].set_title('Error Trend by Total Bit Budget')
        for config in df['config'].unique():
            config_data = df[df['config'] == config].sort_values('total_bits')
            # Group by total_bits and calculate mean error
            trend_data = config_data.groupby('total_bits')['error_value'].agg(['mean', 'std']).reset_index()
            axes[0, 0].plot(trend_data['total_bits'], trend_data['mean'], 
                           marker='o', label=config, linewidth=2, markersize=6)
            # Add error bars
            axes[0, 0].fill_between(trend_data['total_bits'], 
                                   trend_data['mean'] - trend_data['std'],
                                   trend_data['mean'] + trend_data['std'],
                                   alpha=0.2)
        axes[0, 0].set_xlabel('Total Bit Budget')
        axes[0, 0].set_ylabel('Mean Error Value')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # 2. Error trend by channel bit allocation patterns
        axes[0, 1].set_title('Error Trend by Channel Bit Allocation')
        # Create a mapping from channel_bits to a numeric order based on variance
        channel_patterns = sorted(df['channel_bits'].unique(), key=lambda x: np.var([int(x[0]), int(x[1]), int(x[2])]))
        pattern_to_num = {pattern: i for i, pattern in enumerate(channel_patterns)}
        
        for config in df['config'].unique():
            config_data = df[df['config'] == config]
            # Group by channel_bits and calculate mean error
            trend_data = config_data.groupby('channel_bits')['error_value'].mean().reset_index()
            trend_data['pattern_num'] = trend_data['channel_bits'].map(pattern_to_num)
            trend_data = trend_data.sort_values('pattern_num')
            
            axes[0, 1].plot(trend_data['pattern_num'], trend_data['error_value'], 
                           marker='s', label=config, linewidth=2, markersize=6)
        
        axes[0, 1].set_xlabel('Channel Bit Allocation Pattern')
        axes[0, 1].set_ylabel('Mean Error Value')
        axes[0, 1].set_xticks(range(len(channel_patterns)))
        axes[0, 1].set_xticklabels(channel_patterns, rotation=45)
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. Error trend comparison: RGB vs YUV
        axes[1, 0].set_title('Error Trend: RGB vs YUV Color Spaces')
        for colorspace in df['colorspace'].unique():
            colorspace_data = df[df['colorspace'] == colorspace].sort_values('total_bits')
            trend_data = colorspace_data.groupby('total_bits')['error_value'].agg(['mean', 'std']).reset_index()
            axes[1, 0].plot(trend_data['total_bits'], trend_data['mean'], 
                           marker='^', label=colorspace, linewidth=2, markersize=6)
            axes[1, 0].fill_between(trend_data['total_bits'], 
                                   trend_data['mean'] - trend_data['std'],
                                   trend_data['mean'] + trend_data['std'],
                                   alpha=0.2)
        axes[1, 0].set_xlabel('Total Bit Budget')
        axes[1, 0].set_ylabel('Mean Error Value')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # 4. Error trend comparison: equal vs nonuniform quantization
        axes[1, 1].set_title('Error Trend: Equal vs Nonuniform Quantization')
        for quantize_method in df['quantize_method'].unique():
            quantize_data = df[df['quantize_method'] == quantize_method].sort_values('total_bits')
            trend_data = quantize_data.groupby('total_bits')['error_value'].agg(['mean', 'std']).reset_index()
            axes[1, 1].plot(trend_data['total_bits'], trend_data['mean'], 
                           marker='d', label=quantize_method, linewidth=2, markersize=6)
            axes[1, 1].fill_between(trend_data['total_bits'], 
                                   trend_data['mean'] - trend_data['std'],
                                   trend_data['mean'] + trend_data['std'],
                                   alpha=0.2)
        axes[1, 1].set_xlabel('Total Bit Budget')
        axes[1, 1].set_ylabel('Mean Error Value')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('error_trend_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("Trend charts saved as: error_trend_analysis.png")
    
    def create_channel_allocation_trends(self, df):
        """Create line charts showing error trends for same N (total bit budget) with different channel allocations"""
        print("\nGenerating channel allocation trend charts...")
        
        if df.empty:
            print("Error: DataFrame is empty, cannot generate channel allocation charts")
            return
        
        # Set chart style
        plt.style.use('default')
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Error Trends by Channel Bit Allocation (Same Total Bit Budget)', fontsize=16, fontweight='bold')
        
        # Get unique total bit budgets
        total_bits_list = sorted(df['total_bits'].unique())
        
        # 1. N=4 channel allocation trends
        if 4 in total_bits_list:
            n4_data = df[df['total_bits'] == 4].sort_values('channel_bits')
            for config in n4_data['config'].unique():
                config_data = n4_data[n4_data['config'] == config]
                # Create numeric x-axis for channel patterns sorted by variance
                channel_patterns = sorted(config_data['channel_bits'].unique(), 
                                        key=lambda x: np.var([int(x[0]), int(x[1]), int(x[2])]))
                pattern_to_num = {pattern: i for i, pattern in enumerate(channel_patterns)}
                config_data = config_data.copy()
                config_data['pattern_num'] = config_data['channel_bits'].map(pattern_to_num)
                config_data = config_data.sort_values('pattern_num')
                
                axes[0, 0].plot(config_data['pattern_num'], config_data['error_value'], 
                               marker='o', label=config, linewidth=2, markersize=6)
            
            axes[0, 0].set_title('Error Trend for N=4 (Total Bit Budget = 4)')
            axes[0, 0].set_xlabel('Channel Bit Allocation Pattern')
            axes[0, 0].set_ylabel('Error Value')
            axes[0, 0].set_xticks(range(len(channel_patterns)))
            axes[0, 0].set_xticklabels(channel_patterns, rotation=45)
            axes[0, 0].legend()
            axes[0, 0].grid(True, alpha=0.3)
        
        # 2. N=6 channel allocation trends
        if 6 in total_bits_list:
            n6_data = df[df['total_bits'] == 6].sort_values('channel_bits')
            for config in n6_data['config'].unique():
                config_data = n6_data[n6_data['config'] == config]
                channel_patterns = sorted(config_data['channel_bits'].unique(), 
                                        key=lambda x: np.var([int(x[0]), int(x[1]), int(x[2])]))
                pattern_to_num = {pattern: i for i, pattern in enumerate(channel_patterns)}
                config_data = config_data.copy()
                config_data['pattern_num'] = config_data['channel_bits'].map(pattern_to_num)
                config_data = config_data.sort_values('pattern_num')
                
                axes[0, 1].plot(config_data['pattern_num'], config_data['error_value'], 
                               marker='s', label=config, linewidth=2, markersize=6)
            
            axes[0, 1].set_title('Error Trend for N=6 (Total Bit Budget = 6)')
            axes[0, 1].set_xlabel('Channel Bit Allocation Pattern')
            axes[0, 1].set_ylabel('Error Value')
            axes[0, 1].set_xticks(range(len(channel_patterns)))
            axes[0, 1].set_xticklabels(channel_patterns, rotation=45)
            axes[0, 1].legend()
            axes[0, 1].grid(True, alpha=0.3)
        
        # 3. N=8 channel allocation trends
        if 8 in total_bits_list:
            n8_data = df[df['total_bits'] == 8].sort_values('channel_bits')
            for config in n8_data['config'].unique():
                config_data = n8_data[n8_data['config'] == config]
                channel_patterns = sorted(config_data['channel_bits'].unique(), 
                                        key=lambda x: np.var([int(x[0]), int(x[1]), int(x[2])]))
                pattern_to_num = {pattern: i for i, pattern in enumerate(channel_patterns)}
                config_data = config_data.copy()
                config_data['pattern_num'] = config_data['channel_bits'].map(pattern_to_num)
                config_data = config_data.sort_values('pattern_num')
                
                axes[1, 0].plot(config_data['pattern_num'], config_data['error_value'], 
                               marker='^', label=config, linewidth=2, markersize=6)
            
            axes[1, 0].set_title('Error Trend for N=8 (Total Bit Budget = 8)')
            axes[1, 0].set_xlabel('Channel Bit Allocation Pattern')
            axes[1, 0].set_ylabel('Error Value')
            axes[1, 0].set_xticks(range(len(channel_patterns)))
            axes[1, 0].set_xticklabels(channel_patterns, rotation=45)
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)
        
        # 4. Combined view: All N values with channel allocation patterns
        axes[1, 1].set_title('Combined View: All N Values')
        colors = ['red', 'blue', 'green', 'orange']
        markers = ['o', 's', '^', 'd']
        
        for i, n_value in enumerate(total_bits_list):
            n_data = df[df['total_bits'] == n_value]
            # Calculate mean error for each channel pattern across all configs
            mean_errors = n_data.groupby('channel_bits')['error_value'].mean().reset_index()
            channel_patterns = sorted(mean_errors['channel_bits'].unique(), 
                                    key=lambda x: np.var([int(x[0]), int(x[1]), int(x[2])]))
            pattern_to_num = {pattern: j for j, pattern in enumerate(channel_patterns)}
            mean_errors['pattern_num'] = mean_errors['channel_bits'].map(pattern_to_num)
            mean_errors = mean_errors.sort_values('pattern_num')
            
            axes[1, 1].plot(mean_errors['pattern_num'], mean_errors['error_value'], 
                           marker=markers[i % len(markers)], 
                           color=colors[i % len(colors)],
                           label=f'N={n_value}', linewidth=2, markersize=6)
        
        axes[1, 1].set_xlabel('Channel Bit Allocation Pattern')
        axes[1, 1].set_ylabel('Mean Error Value')
        # Use the most comprehensive set of channel patterns for x-axis, sorted by variance
        all_patterns = sorted(df['channel_bits'].unique(), 
                            key=lambda x: np.var([int(x[0]), int(x[1]), int(x[2])]))
        axes[1, 1].set_xticks(range(len(all_patterns)))
        axes[1, 1].set_xticklabels(all_patterns, rotation=45)
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('channel_allocation_trends.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("Channel allocation trend charts saved as: channel_allocation_trends.png")
    
    def generate_detailed_report(self, df):
        """Generate detailed analysis report"""
        print("\n" + "="*60)
        print("Detailed Analysis Report")
        print("="*60)
        
        if df.empty:
            print("Error: DataFrame is empty, cannot generate report")
            return None, None, None, None, None, None
        
        # 1. Configuration comparison
        print("\n1. Configuration Comparison:")
        print("-" * 40)
        config_stats = df.groupby('config')['error_value'].agg(['count', 'mean', 'std', 'min', 'max'])
        config_stats.columns = ['Count', 'Mean Error', 'Std Dev', 'Min Error', 'Max Error']
        print(config_stats.round(0))
        
        # 2. Best configuration analysis
        print("\n2. Best Configuration Analysis:")
        print("-" * 40)
        best_config = df.loc[df['error_value'].idxmin()]
        print(f"Lowest Error Value: {best_config['error_value']:,}")
        print(f"Best Configuration: {best_config['config']}")
        print(f"Color Space: {best_config['colorspace']}")
        print(f"Quantization Method: {best_config['quantize_method']}")
        print(f"Total Bit Budget: {best_config['total_bits']}")
        print(f"Channel Bit Allocation: {best_config['channel_bits']}")
        
        # 3. Color space analysis
        print("\n3. Color Space Analysis:")
        print("-" * 40)
        colorspace_stats = df.groupby('colorspace')['error_value'].agg(['mean', 'std'])
        colorspace_stats.columns = ['Mean Error', 'Std Dev']
        print(colorspace_stats.round(0))
        
        # 4. Quantization method analysis
        print("\n4. Quantization Method Analysis:")
        print("-" * 40)
        quantize_stats = df.groupby('quantize_method')['error_value'].agg(['mean', 'std'])
        quantize_stats.columns = ['Mean Error', 'Std Dev']
        print(quantize_stats.round(0))
        
        # 5. Bit budget analysis
        print("\n5. Bit Budget Analysis:")
        print("-" * 40)
        bit_stats = df.groupby('total_bits')['error_value'].agg(['mean', 'std'])
        bit_stats.columns = ['Mean Error', 'Std Dev']
        print(bit_stats.round(0))
        
        # 6. Channel bit allocation analysis
        print("\n6. Channel Bit Allocation Analysis (Top 10):")
        print("-" * 40)
        channel_stats = df.groupby('channel_bits')['error_value'].agg(['mean', 'count']).sort_values('mean')
        channel_stats.columns = ['Mean Error', 'Count']
        print(channel_stats.head(10).round(0))
        
        # 7. Correlation analysis
        print("\n7. Correlation Analysis:")
        print("-" * 40)
        correlation_data = df[['total_bits', 'q1', 'q2', 'q3', 'error_value']].corr()
        print("Correlation with Error Value:")
        print(correlation_data['error_value'].sort_values(ascending=False))
        
        return config_stats, best_config, colorspace_stats, quantize_stats, bit_stats, channel_stats
    
    def run_analysis(self):
        """Run complete analysis"""
        print("Starting Image Error Analysis...")
        print("="*60)
        
        # Collect data
        self.collect_data()
        
        # Analyze trends
        df = self.analyze_error_trends()
        
        # Create visualizations
        self.create_visualizations(df)
        
        # Create trend line charts
        self.create_trend_visualizations(df)
        
        # Create channel allocation trend charts
        self.create_channel_allocation_trends(df)
        
        # Generate detailed report
        self.generate_detailed_report(df)
        
        print("\n" + "="*60)
        print("Analysis Complete!")
        print("="*60)
        
        return df

def main():
    """Main function"""
    analyzer = ImageErrorAnalyzer()
    df = analyzer.run_analysis()
    
    # Save data to CSV
    if not df.empty:
        df.to_csv('image_error_analysis_data.csv', index=False, encoding='utf-8-sig')
        print("Data saved to: image_error_analysis_data.csv")

if __name__ == "__main__":
    main()