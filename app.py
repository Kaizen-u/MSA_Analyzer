import os
import subprocess
import tempfile
import io
import pandas as pd
import streamlit as st
from Bio import AlignIO, SeqIO


# --- CORE BIOINFORMATICS FUNCTIONS ---

def run_muscle_msa(input_path, output_path):
    """Runs the system-level MUSCLE command."""
    # NOTE: Using -in and -out for MUSCLE v3.8 (Default on Streamlit Cloud's Debian server)
    # If testing locally on Windows with MUSCLE v5, change these to "-align" and "-output"
    cmd = ["muscle", "-in", input_path, "-out", output_path]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        st.error(f"MUSCLE Error: {e.stderr.decode()}")
        st.stop()
    except FileNotFoundError:
        st.error("MUSCLE is not installed. Ensure packages.txt contains 'muscle'.")
        st.stop()


def extract_mutations_horizontal(input_fasta_path, msa_file_path, wt_id):
    """Creates the horizontal mutation matrix, preserving input order."""
    # 1. Grab the exact order of IDs from the ORIGINAL unaligned file
    original_ids = [record.id for record in SeqIO.parse(input_fasta_path, "fasta")]

    # 2. Read the newly aligned file
    alignment = AlignIO.read(msa_file_path, "fasta")

    wt_record, wt_actual_id = None, None
    for record in alignment:
        if record.id == wt_id or wt_id in record.id:
            wt_record = record
            wt_actual_id = record.id
            break

    if not wt_record:
        st.error(f"Wildtype ID '{wt_id}' was not found in the uploaded file.")
        st.stop()

    align_len = alignment.get_alignment_length()
    columns = ["Sequence"] + list(range(1, align_len + 1))
    rows_dict = {}

    # 3. Calculate mutations for all sequences
    for record in alignment:
        row = [record.id]
        if record.id == wt_actual_id:
            row.extend(list(wt_record.seq))
        else:
            for i in range(align_len):
                wt_char = wt_record.seq[i]
                mut_char = record.seq[i]
                # Leave blank if it matches the wildtype
                row.append("" if mut_char == wt_char else mut_char)
        rows_dict[record.id] = row

    # 4. Reorder the final data to match the original input file
    data = []
    for orig_id in original_ids:
        if orig_id in rows_dict:
            data.append(rows_dict[orig_id])
        else:
            # Fallback for slight mismatches caused by Biopython ID chopping
            for msa_id in rows_dict.keys():
                if msa_id in orig_id or orig_id in msa_id:
                    data.append(rows_dict[msa_id])
                    break

    return pd.DataFrame(data, columns=columns)


# --- STREAMLIT UI ---

st.set_page_config(page_title="MSA Mutation Analyzer", layout="wide")
st.title("🧬 MSA Mutation Matrix Generator")
st.markdown("Upload a FASTA file, specify your Wildtype, and generate a horizontal mutation matrix.")

# 1. Upload the File
uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "fas"])

if uploaded_file is not None:
    # Quick parse just to get headers for the dropdown menu
    stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    headers = [record.id for record in SeqIO.parse(stringio, "fasta")]

    # 2. Select the Wildtype
    st.write("---")
    wt_selection = st.selectbox("Select the Wildtype (Reference) Sequence:", headers)

    # 3. The Run Button
    if st.button("Generate Mutation Matrix"):
        with st.spinner("Running MUSCLE alignment and extracting mutations..."):

            # Create hidden temp files for MUSCLE to use
            with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp_in, \
                    tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp_out:

                # Write the uploaded file from RAM to the temp file
                tmp_in.write(uploaded_file.getvalue())
                input_path = tmp_in.name
                output_path = tmp_out.name

            try:
                # Run the logic
                run_muscle_msa(input_path, output_path)
                df_matrix = extract_mutations_horizontal(input_path, output_path, wt_selection)

                # Show a preview in the app
                st.success("Analysis Complete!")
                st.write("### Matrix Preview")
                st.dataframe(df_matrix.head(10))  # Show first 10 rows

                # Convert DataFrame to an Excel file in RAM
                output_excel = io.BytesIO()
                with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
                    df_matrix.to_excel(writer, index=False, sheet_name='Mutations')
                processed_data = output_excel.getvalue()

                # Provide the Download Button
                st.download_button(
                    label="📥 Download Excel Matrix",
                    data=processed_data,
                    file_name="MUT_Ana_horizontal.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )

            finally:
                # Clean up the temp files so the server's hard drive doesn't fill up
                if os.path.exists(input_path): os.remove(input_path)
                if os.path.exists(output_path): os.remove(output_path)