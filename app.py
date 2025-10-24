import streamlit as st
from streamlit_option_menu import option_menu
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from PIL import Image
import io
import base64
import glob
import os

# Page configuration
st.set_page_config(
    page_title="DNABERT-Enhancer Portal",
    page_icon="🧬",
    layout="wide"
)

# Custom CSS
st.markdown("""
<style>
    .main {
        background-color: #f8f9fa;
    }
    .stApp {
        background-color: #f8f9fa;
    }
    div[data-testid="stMetricValue"] {
        font-size: 2rem;
        font-weight: 600;
    }
    h1 {
        color: #1a1a1a;
        font-weight: 700;
    }
    .subtitle {
        color: #6b7280;
        font-size: 1.1rem;
        margin-bottom: 2rem;
    }
    .info-box {
        background: #e0f2fe;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #0284c7;
        margin: 1rem 0;
    }
    .clinical-tag {
        background: #dcfce7;
        color: #166534;
        padding: 0.25rem 0.5rem;
        border-radius: 4px;
        font-size: 0.85rem;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)

# Load data function
@st.cache_data
def load_data():
    detail_datasets = {
        'Enhancer GOF': pd.read_csv('./data/Detailed_info_enhancer_GOF.csv'),
        'Enhancer LOF': pd.read_csv('./data/Detailed_info_enhancer_LOF.csv'),
        'Non-enhancer GOF': pd.read_csv('./data/Detailed_info_non-enhancer_GOF.csv')
    }
    combined_df = pd.concat(detail_datasets.values(), ignore_index=True)
    return combined_df
    

# Sidebar
with st.sidebar:
    st.markdown("## 🧬 DNABERT-Enhancer Portal")
    st.markdown("**An interactive platform for exploration of predictions by DNABERT-Enhancer-350 model**")
    st.divider()

    with st.sidebar:
        st.markdown("### Navigation")
        page = st.radio(
            "",
            ["ℹ️ About", "📊 Browse Data"],
            index=0,
            label_visibility="collapsed",
            horizontal=False
        )
  
# Main content
if page == "📊 Browse Data":
    # Header
    st.title("DNABERT-Enhancer prediction data")

    tab1, tab2 = st.tabs(["Candidate Variants", "Enhancers in Human Genome"])
    with tab1:
        st.markdown(
            """
            <h1 style="font-size:20px; font-weight:bold; color:#1f2937;">
                Candidate Variants predicted by DNABERT-Enhancer
            </h1>
            """,
            unsafe_allow_html=True
        )
        
        combined_df = load_data()
        # st.write("Columns in combined_df:", combined_df.columns.tolist())
        columns_order = [
            "ID","chromosome","region_coordinates","dbsnp_id","variant_start","variant_end","reference_nucleotide",
            "alternative_nucleotide","reference_probability","alternative_probability","ScoreChange","LogOddRatio",
            "reported_clinical_association","gwas_url","clinvar_url","eqtl_url","predicted_functional_effect","class",
            "transcription_factor","tf_reference_probability","tf_alternative_probability","tf_ScoreChange","tf_LogOddRatio",
            "gene","strand","distance","element_coordinates","variant_coordinates"
        ]
        combined_df = combined_df[columns_order]

        col1, col2 = st.columns([1, 2], gap="large")

        with col1:
            st.markdown(
            """
            <h1 style="font-size:18px; font-weight:bold; color:#1f2937;">
                Filter data
            </h1>
            """,
            unsafe_allow_html=True
        )
            
            # Initialize persistent states
            if "filter_key" not in st.session_state:
                st.session_state.filter_key = 0

            if "lor_range" not in st.session_state:
                st.session_state.lor_range = (
                    float(combined_df["LogOddRatio"].min()),
                    float(combined_df["LogOddRatio"].max())
                )

            # Dropdown options
            effect_options = ["All"] + sorted(combined_df["predicted_functional_effect"].dropna().unique())
            class_options = ["All"] + sorted(combined_df["class"].dropna().unique())
            chrom_options = ["All"] + sorted(combined_df["chromosome"].dropna().unique())
            assoc_options = ["All"] + sorted(combined_df["reported_clinical_association"].dropna().unique())

            # --- Widgets (no direct assignment to session_state) ---
            selected_effect = st.selectbox(
                "Predicted Functional Effect",
                effect_options,
                key=f"predicted_effect_{st.session_state.filter_key}"
            )

            selected_class = st.selectbox(
                "Variant Class",
                class_options,
                key=f"variant_class_{st.session_state.filter_key}"
            )

            st.markdown("---")
            st.markdown(
            """
            <h1 style="font-size:18px; font-weight:bold; color:#1f2937;">
                Additional filter
            </h1>
            """,
            unsafe_allow_html=True
        )

            
            selected_chrom = st.selectbox(
                "Chromosome",
                chrom_options,
                key=f"chrom_{st.session_state.filter_key}"
            )

            selected_assoc = st.selectbox(
                "Clinical Association",
                assoc_options,
                key=f"assoc_{st.session_state.filter_key}"
            )

            # Slider
            min_lor, max_lor = float(combined_df["LogOddRatio"].min()), float(combined_df["LogOddRatio"].max())
            lor_range = st.slider(
                "LogOddRatio range",
                min_lor, max_lor,
                value=st.session_state.lor_range,
                key=f"lor_slider_{st.session_state.filter_key}"
            )

            # Reset filters
            if st.button("🔄 Reset Filters"):
                st.session_state.filter_key += 1
                st.session_state.lor_range = (float(min_lor), float(max_lor))
                st.rerun()

            # --- FILTERING LOGIC ---
            filtered_df = combined_df.copy()

            if selected_effect != "All":
                filtered_df = filtered_df[filtered_df["predicted_functional_effect"] == selected_effect]
            if selected_class != "All":
                filtered_df = filtered_df[filtered_df["class"] == selected_class]
            if selected_chrom != "All":
                filtered_df = filtered_df[filtered_df["chromosome"] == selected_chrom]
            if selected_assoc != "All":
                filtered_df = filtered_df[filtered_df["reported_clinical_association"] == selected_assoc]

            min_lor_val, max_lor_val = lor_range
            filtered_df = filtered_df[
                (filtered_df["LogOddRatio"] >= min_lor_val) &
                (filtered_df["LogOddRatio"] <= max_lor_val)
            ]

            # unique_key_cols = ["chromosome", "dbsnp_id", "variant_start", "variant_end", "reference_nucleotide", "alternative_nucleotide"]

            # # Count unique rows based on these columns
            # num_unique_variants = len(filtered_df.drop_duplicates(subset=unique_key_cols))
            # st.markdown(f"**{num_unique_variants:,} variants displayed**")

        with col2:
    		# --- Search Bar and Clear Button ---
            search_col, clear_col = st.columns([5, 0.6])

            if "search_query" not in st.session_state:
                st.session_state.search_query = ""

            with search_col:
                search_query = st.text_input(
                    "🔍 Search Variants",
                    placeholder="Search by rsID, CV ID, chromosome, or keyword...",
                    value=st.session_state.search_query,
                    key=f"search_{st.session_state.filter_key}"
                )

            with clear_col:
                if st.button("❌", help="Clear Search"):
                    st.session_state.search_query = ""
                    st.session_state.filter_key += 1
                    st.rerun()

            display_cols = [
                "ID", "chromosome", "dbsnp_id", "ScoreChange", "LogOddRatio",
                "reported_clinical_association", "predicted_functional_effect", "class"
            ]

            # --- Apply Search Filter ---
            filtered_display_df = filtered_df.copy()
            if search_query:
                search_query = search_query.lower().strip()
                filtered_display_df = filtered_display_df[
                    filtered_display_df.astype(str).apply(
                        lambda row: row.str.lower().str.contains(search_query, na=False)
                    ).any(axis=1)
                ]

            filtered_display_df = filtered_display_df[display_cols].drop_duplicates()

            st.dataframe(filtered_display_df[display_cols], use_container_width=True, height=500, hide_index=True)

            # --- Download option ---
            csv = filtered_display_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Download Filtered Variants (CSV)",
                data=csv,
                file_name="filtered_candidate_variants.csv",
                mime="text/csv"
            )
        
        st.markdown("---")
        st.markdown(
            """
            <h1 style="font-size:20px; font-weight:bold; color:#1f2937;">
                Detailed information on candidate variants
            </h1>
            """,
            unsafe_allow_html=True
        )
        selected_variant_id = st.selectbox(
                "Select Candidate Variant ID to see details",
                options=filtered_display_df["ID"].unique()
            )
                
    with tab2:
        # Load and combine all split files
        data_path = "./data/whole_genome_prediction_data/"  # change to the folder where your CSVs are
        all_files = sorted(glob.glob(os.path.join(data_path, "WGP_*.csv")))
        # Combine all parts into one DataFrame
        df_list = [pd.read_csv(f) for f in all_files]
        combined_wgp_df = pd.concat(df_list, ignore_index=True)
        
        st.header("Enhancers in Human Genome")
        
        st.dataframe(combined_wgp_df, use_container_width=True, height=500, hide_index=True)
        
else:  # About page
    st.title("About DNABERT-Enhancer portal")
    st.markdown("""
    <div style="text-align: justify;">

    <h5>Overview:</h5>
    The DNABERT-Enhancer portal offers an interactive platform to explore candidate gain- and loss-of-function enhancer variants predicted by 
    the DNABERT-Enhancer-350 model using ENCODE SCREEN enhancers (350 bp). It also provides access to genome-wide enhancer predictions across 
    the human reference genome (GRCh38). This web application enables users to visualize, search, and interpret enhancer regions and their 
    potential functional impact in a genomic context.
    </div>
    """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    
    # Load the image
    img = Image.open("./Figures/Enhancer.png")

    # Save to a buffer to preserve quality
    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    buf.seek(0)

    # Encode image in base64
    base64_img = base64.b64encode(buf.read()).decode("utf-8")

    left_col, right_col = st.columns([2, 1], gap="large")
    with left_col:
        st.markdown("""
        <div style="text-align: justify;">
        <h5>Background:</h5>
        Enhancers are one of the cis-regulatory element which increases the transcription of a target genes while interacting with their target promoters with the assistance of proteins like transcription factors, 
        mediators and RNA polymerase, thereby shaping the characteristics and function of cells and tissues. Any disruption in the ideal function of enhancers due to genetic or epigenetic changes leads to disease conditions. 
        Predicting and deciphering the regulatory logic of enhancers is a challenging problem, due to the intricate sequence features and lack of consistent genetic or epigenetic signatures that can accurately discriminate enhancers 
        from other genomic regions. <b>DNABERT-Enhancer</b>, a novel enhancer prediction method developed by applying DNABERT pre-trained language model on the human genome, learns the "language of DNA" and predict enhancer activity directly 
        from sequence, providing genome-wide enhancer annotation and variant impact prediction. These genome-wide enhancers and candidate genetic variants predicted by DNABERT-Enhancer provide valuable resources for 
        genome interpretation in functional and clinical genomics studies.
        </div>
        """, unsafe_allow_html=True)

    with right_col:
        # Display image with padding from top
        st.markdown(
        f"""
        <div style="padding-top: 40px; text-align: center;">
            <h6 style="margin-bottom: 15px;">
                Enhancer–promoter interaction mediated by transcription factors, mediator complex, and RNA polymerase II
            </h6>
            <img src="data:image/png;base64,{base64_img}" width="300">
        </div>
        """,
        unsafe_allow_html=True
        )

    img2 = Image.open("./Figures/Graphical_abstract.png")
    # Save to buffer to preserve quality
    buf = io.BytesIO()
    img2.save(buf, format="PNG", optimize=True)
    buf.seek(0)
    # Encode to base64
    base64_img2 = base64.b64encode(buf.read()).decode("utf-8")
    # Display centered with small width
    st.markdown(
        f"""
        <div style="text-align: center; padding-top: 10px;">
            <img src="data:image/png;base64,{base64_img2}" width="600">
            <h6 style="margin-bottom: 15px;">
                DNABERT-Enhancer Model : Study outline
            </h6>
        </div>
        """,
        unsafe_allow_html=True
    )
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    st.markdown("""
    <div style="text-align: justify;">
    <h5>Datasets available:</h5>
    <h6>1. Candidate variants:</h6>
    The DNABERT-Enhancer identifies candidate regulatory variants by evaluating short variants from dbSNP release 155 (GRCh38) located 
    within enhancer regions and transcription factor (TF) target sites. Variant effects are predicted by substituting alternate alleles into 
    enhancer or TF-target sites within and re-evaluating them using DNABERT-Enhancer-350 and DeepVRegulome, respectively. Variants are classified 
    as gain-of-function (GOF) or loss-of-function (LOF) based on changes in model prediction probabilities, score differences, and log-odds 
    ratios (LOR). Significance thresholds are defined empirically from the LOR distribution to highlight variants with strong regulatory effects. 
    The data also includes functional annotations from ClinVar, GWAS Catalog, and GTEx eQTL data to link predicted effects with known 
    clinical and expression traits along with information on nearest gene.
    </div>
        """, unsafe_allow_html=True)
    
    data = {
        "Dataset": ["Enhancer GOF", "Enhancer LOF", "Non-enhancer GOF"],
        "Number of Variants": [1917, 2681, 5464],
        "Description": [
            "Variants that increase enhancer activity within known enhancers",
            "Variants that decrease enhancer activity within known enhancers; affecting 301 TFs",
            "Variants that create new enhancer activity outside canonical enhancers",
        ],
        "Score Interpretation": [
            "Positive score changes indicate gain of function",
            "Negative score changes indicate loss of function",
            "Positive score changes indicate gain of function",
        ]
    }
    df = pd.DataFrame(data)
    st.markdown(
        df.to_html(index=False, justify="center"),
        unsafe_allow_html=True
    )

    st.markdown("""
    <div style="text-align: justify;">
    <h6>2. Enhancers in Human Genome:</h6>
    The DNABERT-Enhancer-350 model is applied to the human reference genome (GRCh38) to generate genome-wide enhancer predictions. 
    The genome is segmented into 350 bp sequences using a 150 bp sliding window, excluding regions with unidentified bases (‘N’). 
    Each sequence is tokenized into hexamers and evaluated by the model, labeling windows with prediction probabilities ≥ 0.5 as 
    enhancers. Consecutive enhancer-positive windows are merged into continuous regions, resulting in the identification of 1.8M 
    enhancer regions, collectively covering 21.53% of the human genome. The predicted enhancers are compared with known enhancer 
    databases at multiple overlap thresholds (5–95%) to assess concordance, and statistical significance is evaluated through permutation 
    testing.
    </div>
        """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    
    st.markdown("""
    <div style="text-align: justify;">
        <h5>Features:</h5>
        <ul>
            <li>Advanced filtering by chromosome, position, and score change</li>
            <li>Clinical significance filtering</li>
            <li>External database links (GWAS, ClinVar, eQTL)</li>
            <li>Transcription factor analysis</li>
            <li>Comparative visualizations</li>
            <li>Downloadable filtered results</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    st.markdown("""
    <div style="text-align: justify;">
        <h5>Interpretation:</h5>
        <ul>
            <li><b>ScoreChange:</b> Magnitude of functional impact</li>
            <li><b>LogOddRatio:</b> Statistical confidence in prediction</li>
            <li><b>Clinical Significance:</b> Known disease associations</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    
    st.markdown("""
    <div style="text-align: justify;">
    <h5>Intended use:</h5>
    This portal is designed for research and exploratory use only. Predictions are computational and should be interpreted as hypothesis-generating, 
    not diagnostic or clinical evidence.
    </div>
        """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    
    st.markdown("""
    <div style="text-align: justify;">
    <h5>Citation:</h5>
    If you use DNABERT-Enhancer or data from this portal, please cite:<br>
    Sathian R, Dutta P, Ay F, Davuluri RV. <b>Genomic Language Model for Predicting Enhancers and Their Allele-Specific Activity in the Human Genome</b>, 2025.  
    <a href="https://doi.org/10.1101/2025.03.18.644040" target="_blank">DOI: 10.1101/2025.03.18.644040</a>
    </div>
        """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    
    st.markdown("""
    <div style="text-align: justify;">
    <h5>Acknowledgements:</h5>
    We thank all members of the Davuluri lab ( 548 The State University of New York at Stony Brook) and Dante Bolzan 
    (Ay Lab - La Jolla Institute for Immunology) for critical discussions and helpful 
    advice. This work was financially supported by grants from National Library of Medicine/National
    Institutes of Health funding – [R01LM01372201 to R.D., R35GM128938 to F.A]<br>
    For inquiries or collaborations, please contact: Tel: +1 (631) 638-2590; Email: Ramana.Davuluri@stonybrookmedicine.edu.<br>
    <br>
    Built with Streamlit 🎈
    </div>
        """, unsafe_allow_html=True)
