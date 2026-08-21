#!/usr/bin/env python3

import pandas as pd
import json
import sys
import os

def load_clinical_metadata(file_path):

    if not os.path.exists(file_path):
        print(f"Clinical metadata file not found: {file_path}")
        return {}
    
    try:
        if file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path)
        elif file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        else:
            print(f"Unsupported file format: {file_path}")
            return {}
        
        metadata_dict = {}
        for _, row in df.iterrows():
            sample_id = str(row.get('sample_id', row.get('Sample_ID', row.get('SampleID', ''))))
            if sample_id:
                metadata_dict[sample_id] = row.to_dict()
        
        return metadata_dict
    
    except Exception as e:
        return {}

def find_matching_sample(sample_id, metadata_dict):

    if sample_id in metadata_dict:
        return metadata_dict[sample_id]
    
    variations = [
        sample_id.replace('_ont', ''),
        sample_id.replace('_illumina', ''),
        sample_id.replace('.fastq', ''),
        sample_id.replace('.consensus', ''),
        sample_id.replace('_1', '').replace('_2', ''),
        sample_id.split('_')[0] if '_' in sample_id else sample_id
    ]

    for variation in variations:
        if variation in metadata_dict:
            return metadata_dict[variation]

    return None

def get_clinical_value(metadata, field, default='Unknown'):
    if not metadata:
        return default
    
    possible_fields = [
        field,
        field.lower(),
        field.upper(),
        field.replace('_', ' ').title(),
        field.replace(' ', '_').lower()
    ]
    
    for f in possible_fields:
        if f in metadata and pd.notna(metadata[f]) and str(metadata[f]).strip():
            return str(metadata[f]).strip()
    
    return default


def load_organization_metadata(file_path):
    if not os.path.exists(file_path):
        print(f"Organization metadata file not found: {file_path}")
        return {}
    try:
        if file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path)
        elif file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        else:
            print(f"Unsupported file format: {file_path}")
            return {}
        if df.empty:
            return {}
        row = df.iloc[0].to_dict()
        return {k: (str(v).strip() if pd.notna(v) else '') for k, v in row.items()}
    except Exception as e:
        print(f"Error loading organization metadata: {e}")
        return {}


def load_practitioner_metadata(file_path):
    if not os.path.exists(file_path):
        print(f"Practitioner metadata file not found: {file_path}")
        return {}
    try:
        if file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path)
        elif file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        else:
            print(f"Unsupported file format: {file_path}")
            return {}
        if df.empty:
            return {}
        row = df.iloc[0].to_dict()
        return {k: (str(v).strip() if pd.notna(v) else '') for k, v in row.items()}
    except Exception as e:
        print(f"Error loading practitioner metadata: {e}")
        return {}