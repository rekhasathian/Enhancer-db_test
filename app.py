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

st.markdown("""
<style>
.stTabs [data-testid="stWaiver"] button {
        font-size: 16px;
        font-weight: 700;
        color: #3182CE;
        font-family: Arial, sans-serif;
    }
    .stTabs [data-testid="stWaiver"] button[aria-selected="true"] {
        color: #085D9E;
        border-bottom: 2px solid #085D9E;
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

# If there's a variant in URL, default to Browse
query_params = st.experimental_get_query_params()
variant_in_url = "variant" in query_params and query_params.get("variant")

# Sidebar
with st.sidebar:
    st.markdown("## 🧬 DNABERT-Enhancer Portal")
    st.markdown("**An interactive platform for exploration of predictions by DNABERT-Enhancer-350 model**")
    st.divider()

    st.markdown("### Navigation")
    # default_index selects Browse if variant present
    default_index = 1 if variant_in_url else 0
    page = st.radio(
        "",
        ["ℹ️ About", "📊 Browse Data"],
        index=default_index,
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
        columns_order = [
            "ID","chromosome","region_coordinates","dbsnp_id","variant_start","variant_end","reference_nucleotide",
            "alternative_nucleotide","reference_probability","alternative_probability","ScoreChange","LogOddRatio",
            "reported_clinical_association","gwas_url","clinvar_url","eqtl_url","predicted_functional_effect","class",
            "transcription_factor","tf_reference_probability","tf_alternative_probability","tf_ScoreChange","tf_LogOddRatio",
            "gene","strand","distance","element_coordinates","variant_coordinates"
        ]
        # ensure columns exist before reordering
        combined_df = combined_df[[c for c in columns_order if c in combined_df.columns]]
        col1, col2 = st.columns([1, 4], gap="medium")

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
                if "LogOddRatio" in combined_df.columns:
                    st.session_state.lor_range = (
                        float(combined_df["LogOddRatio"].min()),
                        float(combined_df["LogOddRatio"].max())
                    )
                else:
                    st.session_state.lor_range = (0.0, 1.0)

            # Dropdown options (guard missing columns)
            effect_options = ["All"]
            if "predicted_functional_effect" in combined_df.columns:
                effect_options += sorted(combined_df["predicted_functional_effect"].dropna().unique().tolist())

            class_options = ["All"]
            if "class" in combined_df.columns:
                class_options += sorted(combined_df["class"].dropna().unique().tolist())

            chrom_options = ["All"]
            if "chromosome" in combined_df.columns:
                chrom_options += sorted(combined_df["chromosome"].dropna().unique().tolist())

            assoc_options = ["All"]
            if "reported_clinical_association" in combined_df.columns:
                assoc_options += sorted(combined_df["reported_clinical_association"].dropna().unique().tolist())

            # --- Widgets (no direct assignment to session_state) ---
            selected_effect = st.selectbox(
                "Predicted Functional Effect",
                effect_options,
                key=f"predicted_effect_{st.session_state.filter_key}"
            )

            selected_class = st.selectbox(
                "Element Class",
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
                "Reported Clinical Association",
                assoc_options,
                key=f"assoc_{st.session_state.filter_key}"
            )

            # Slider
            if "LogOddRatio" in combined_df.columns:
                min_lor, max_lor = float(combined_df["LogOddRatio"].min()), float(combined_df["LogOddRatio"].max())
            else:
                min_lor, max_lor = 0.0, 1.0

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

            if selected_effect != "All" and "predicted_functional_effect" in filtered_df.columns:
                filtered_df = filtered_df[filtered_df["predicted_functional_effect"] == selected_effect]
            if selected_class != "All" and "class" in filtered_df.columns:
                filtered_df = filtered_df[filtered_df["class"] == selected_class]
            if selected_chrom != "All" and "chromosome" in filtered_df.columns:
                filtered_df = filtered_df[filtered_df["chromosome"] == selected_chrom]
            if selected_assoc != "All" and "reported_clinical_association" in filtered_df.columns:
                filtered_df = filtered_df[filtered_df["reported_clinical_association"] == selected_assoc]

            min_lor_val, max_lor_val = lor_range
            if "LogOddRatio" in filtered_df.columns:
                filtered_df = filtered_df[
                    (filtered_df["LogOddRatio"] >= min_lor_val) &
                    (filtered_df["LogOddRatio"] <= max_lor_val)
                ]

        with col2:
            # --- Search Bar and Clear Button ---
            search_col, clear_col = st.columns([4, 0.5], gap="small")

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
                st.markdown(
                    """
                    <style>
                    div[data-testid="stButton"] button {
                        padding: 0.2rem 0.4rem;
                        font-size: 0.8rem;
                        margin-top: 0.8rem; /* vertically align with search bar */
                    }
                    </style>
                    """,
                    unsafe_allow_html=True,
                )
                if st.button("❌", help="Clear Search"):
                    st.session_state.search_query = ""
                    st.session_state.filter_key += 1
                    st.rerun()

            display_cols = [
                "ID", "chromosome", "dbsnp_id", "ScoreChange", "LogOddRatio",
                "reported_clinical_association", "predicted_functional_effect", "class"
            ]
            # keep only existing columns
            display_cols = [c for c in display_cols if c in filtered_df.columns]

            # --- Apply Search Filter ---
            filtered_display_df = filtered_df.copy()
            if search_query:
                search_query = search_query.lower().strip()
                filtered_display_df = filtered_display_df[
                    filtered_display_df.astype(str).apply(
                        lambda row: row.str.lower().str.contains(search_query, na=False)
                    ).any(axis=1)
                ]

                if filtered_display_df.empty:
                    st.warning(f"No results found for **'{search_query}'**.")
                
            filtered_display_df = filtered_display_df[display_cols].drop_duplicates()

            # --- Make ID clickable (fast HTML) ---
            if "ID" in filtered_display_df.columns:
                filtered_display_df = filtered_display_df.copy()
                filtered_display_df["ID"] = filtered_display_df["ID"].apply(
                    lambda x: f'<a href="?variant={x}" target="_self" style="color:#0073e6; text-decoration:none;">{x}</a>'
                )

            # --- Optimize table display for speed ---
            max_rows = 50  # number of rows to display (you can make it a user option too)

            # Apply all filters/search as before
            # filtered_display_df = <your existing filtered dataframe>

            # Limit display but keep all filter logic intact
            display_subset = filtered_display_df.head(max_rows)

            st.markdown("""
            <style>
            .scroll-table-container {
                max-height: 500px;
                overflow-y: auto;
                overflow-x: auto;
                border-bottom: 2px solid #ddd;
                background-color: #ffffff;
                border-radius: 6px;
            }
            .scroll-table-container table {
                border-collapse: collapse;
                width: 100%;
                font-size: 13px;
            }
            .scroll-table-container thead th {
                position: sticky;
                top: 0;
                background-color: #f3f6fa;
                color: #333;
                font-weight: 600;
                border-bottom: 2px solid #ccc;
                text-align: left;
                z-index: 2;
            }
            .scroll-table-container th, .scroll-table-container td {
                padding: 6px 10px;
                border-bottom: 1px solid #e6e6e6;
            }
            .scroll-table-container tr:nth-child(even) {
            background-color: #fafafa;
            }
            </style>
            """, unsafe_allow_html=True)
            
            st.markdown(
                f"""
                <div class="scroll-table-container">
                    {display_subset.to_html(escape=False, index=False)}
                </div>
                """,
                unsafe_allow_html=True
            )

            # --- Optional: Show row count info ---
            st.caption(f"Showing top {len(display_subset):,} of {len(filtered_display_df):,} matching rows.")


            # small JS: intercept clicks on ?variant links and update query params w/o full reload
            st.markdown(
                """
                <script>
                document.addEventListener('click', function(event) {
                    const link = event.target.closest('a[href^="?variant="]');
                    if (link) {
                        event.preventDefault();  // stop full reload
                        const urlParams = new URLSearchParams(link.getAttribute('href').substring(1));
                        const variant = urlParams.get('variant');
                        // set the query param in Streamlit without a full reload
                        window.parent.postMessage({ type: 'streamlit:setQueryParams', params: { variant } }, '*');
                    }
                });
                </script>
                """,
                unsafe_allow_html=True
            )

            # --- Download option ---
            csv = filtered_display_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Download Filtered Variants (CSV)",
                data=csv,
                file_name="filtered_candidate_variants.csv",
                mime="text/csv"
            )
        
        # --- Detailed info section (below the table) ---
        st.markdown("---")
        
        query_params_now = st.experimental_get_query_params()
        if "variant" in query_params_now:
            selected_variant_id = (
                query_params_now["variant"][0]
                if isinstance(query_params_now["variant"], list)
                else query_params_now["variant"]
            )

            detailed_info = combined_df[combined_df["ID"] == selected_variant_id]
            if not detailed_info.empty:
                row = detailed_info.iloc[0]
                rowd = row.to_dict()

                def pick(*keys):
                    for key in keys:
                        if key in rowd and pd.notna(rowd[key]):
                            return rowd[key]
                    return "N/A"

                st.markdown(f"### 🧬 Detailed information for variant: `{selected_variant_id}`")

                gene_id = pick('gene', 'Closest gene', 'gene_id')
                ensembl_link = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={gene_id}" if gene_id != "N/A" else None
                rs_id = pick('dbsnp_id')
                dbsnp_link = f"https://www.ncbi.nlm.nih.gov/snp/{rs_id}" if rs_id != "N/A" else None

                # --- Custom card layout ---
                st.markdown("""
                <style>
                .info-card {
                    background-color: #ffffff;
                    border-radius: 15px;
                    box-shadow: 0 2px 6px rgba(0,0,0,0.08);
                    padding: 20px;
                    margin-bottom: 20px;
                }
                .info-title {
                    font-size: 18px;
                    font-weight: 600;
                    color: #1f2937;
                    margin-bottom: 10px;
                }
                .info-grid {
                    display: grid;
                    grid-template-columns: 1fr 1fr;
                    gap: 8px 40px;
                }
                .info-label {
                    font-weight: 600;
                    color: #374151;
                }
                .info-value {
                    color: #111827;
                }
                </style>
                """, unsafe_allow_html=True)

                # --- Basic Information Card ---
                st.markdown("""
                <div class="info-card">
                    <div class="info-title"> Basic Information</div>
                    <div class="info-grid">
                        <div class="info-label">Candidate Variant ID</div><div class="info-value">{id}</div>
                        <div class="info-label">Genomic Element Class</div><div class="info-value">{cls}</div>
                        <div class="info-label">Organism</div><div class="info-value">Human</div>
                        <div class="info-label">Genome Assembly</div><div class="info-value">GRCh38</div>
                        <div class="info-label">Element Coordinate</div><div class="info-value">{coord}</div>
                        <div class="info-label">Closest Gene</div><div class="info-value"><a href="{ens}" target="_blank">{gene}</a></div>
                        <div class="info-label">Strand</div><div class="info-value">{strand}</div>
                        <div class="info-label">Distance to Element</div><div class="info-value">{dist}</div>
                    </div>
                </div>
                """.format(
                    id=pick('ID'),
                    cls=pick('class'),
                    coord=pick('element_coordinates'),
                    gene=gene_id,
                    ens=ensembl_link,
                    strand=pick('strand'),
                    dist=pick('distance')
                ), unsafe_allow_html=True)

                # --- Variant Prediction Card ---
                st.markdown("""
                <div class="info-card">
                    <div class="info-title"> Variant Prediction Information</div>
                    <div class="info-grid">
                        <div class="info-label">Reference SNP (rs) ID</div><div class="info-value"><a href="{dbsnp}" target="_blank">{rs}</a></div>
                        <div class="info-label">Variant Coordinate</div><div class="info-value">{varcoord}</div>
                        <div class="info-label">Reference Allele</div><div class="info-value">{ref}</div>
                        <div class="info-label">Alternative Allele</div><div class="info-value">{alt}</div>
                        <div class="info-label">Reference Probability</div><div class="info-value">{refprob}</div>
                        <div class="info-label">Alternative Probability</div><div class="info-value">{altprob}</div>
                        <div class="info-label">Score Change</div><div class="info-value">{sc}</div>
                        <div class="info-label">Log Odds Ratio</div><div class="info-value">{lor}</div>
                    </div>
                </div>
                """.format(
                    dbsnp=dbsnp_link,
                    rs=rs_id,
                    varcoord=pick('variant_coordinates'),
                    ref=pick('reference_nucleotide'),
                    alt=pick('alternative_nucleotide'),
                    refprob=pick('reference_probability'),
                    altprob=pick('alternative_probability'),
                    sc=pick('ScoreChange'),
                    lor=pick('LogOddRatio')
                ), unsafe_allow_html=True)

                st.markdown("##### Reported Clinical Association")

                reported_raw = rowd.get("reported_clinical_association", None)
                reported = str(reported_raw).strip().lower() if reported_raw is not None and not pd.isna(reported_raw) else "no"
                # st.write(f"Reported raw value: '{reported_raw}'")
                
                clinvar_url = rowd.get("clinvar_url", None)
                gwas_url = rowd.get("gwas_url", None)
                eqtl_url = rowd.get("eqtl_url", None)
                
                if reported_raw == "yes":
                    st.markdown("**This variant has reported clinical significance.**")

                    links = []
                    if clinvar_url and str(clinvar_url).strip() not in ["-", "nan", "", "None"]:
                        links.append(f'<a href="{clinvar_url}" target="_blank" style="text-decoration:none;">🧫 <b>ClinVar</b></a>')
                    if gwas_url and str(gwas_url).strip() not in ["-", "nan", "", "None"]:
                        links.append(f'<a href="{gwas_url}" target="_blank" style="text-decoration:none;">📊 <b>GWAS Catalog</b></a>')
                    if eqtl_url and str(eqtl_url).strip() not in ["-", "nan", "", "None"]:
                        links.append(f'<a href="{eqtl_url}" target="_blank" style="text-decoration:none;">🧠 <b>GTEx eQTL</b></a>')

                    if links:
                        links_html = " &nbsp; | &nbsp; ".join(links)

                        st.markdown(
                            f"""
                            <div style="
                                background-color:#f0f8ff;
                                border: 1px solid #b3d4fc;
                                border-radius: 12px;
                                padding: 12px 16px;
                                margin-top: 6px;
                                font-size: 14px;
                            ">
                                🔗 {links_html}
                            </div>
                            """,
                            unsafe_allow_html=True
                        )
                    else:
                        st.markdown(
                            """
                            <div style="
                                background-color:#fff8e1;
                                border: 1px solid #f5d742;
                                border-radius: 12px;
                                padding: 12px 16px;
                                color:#555;
                                font-size:14px;
                            ">
                                ⚠️ Marked as reported, but no external database links provided.
                            </div>
                            """,
                        unsafe_allow_html=True
                    )

                else:
                    st.markdown(
                        """
                        <div style="
                            background-color:#fafafa;
                            border:1px solid #ddd;
                            border-radius:12px;
                            padding:12px 16px;
                            color:#666;
                            font-size:14px;
                        ">
                            ⚪ No reported clinical significance found.
                        </div>
                        """,
                        unsafe_allow_html=True
                    )
                # --- Transcription Factor Binding Impact (for LOF variants only) ---
                pred_effect = str(rowd.get("predicted_functional_effect", "")).strip().lower()
                # st.write(f"predicted functional effect: '{pred_effect}'")
                if "loss of function" in pred_effect:  # show only for LOF variants
                    # Filter all TF rows for this variant
                    tf_rows = combined_df[combined_df["ID"] == selected_variant_id][[
                        "transcription_factor",
                        "tf_reference_probability",
                        "tf_alternative_probability",
                        "tf_ScoreChange",
                        "tf_LogOddRatio"
                    ]].dropna(subset=["transcription_factor"], how="all")

                    if not tf_rows.empty:
                        # Card header
                        st.markdown("##### Transcription Factor Binding impact")

                        # Rename columns for display
                        tf_rows_display = tf_rows.rename(columns={
                            "transcription_factor": "Transcription Factor",
                            "tf_reference_probability": "Reference Probability",
                            "tf_alternative_probability": "Alternative Probability",
                            "tf_ScoreChange": "Score Change",
                            "tf_LogOddRatio": "Log Odds Ratio"
                        })

                        # Optional: clean numeric formatting
                        for col in ["Reference Probability", "Alternative Probability", "Score Change", "Log Odds Ratio"]:
                            if col in tf_rows_display.columns:
                                tf_rows_display[col] = tf_rows_display[col].apply(
                                    lambda x: f"{x:.4f}" if isinstance(x, (float, int)) else x
                                )

                        # Display compact, scrollable table
                        st.dataframe(
                            tf_rows_display.style.set_table_styles([
                                {'selector': 'thead th', 'props': [('background-color', '#f0f8ff'), ('font-weight', '600')]},
                                {'selector': 'tbody td', 'props': [('font-size', '13px')]}
                            ]),
                            use_container_width=True,
                            hide_index=True,
                            height=min(300, 40 + len(tf_rows_display) * 30)
                        )

                    else:
                        st.markdown(
                            """
                            <div style="
                                background-color:#fafafa;
                                border:1px solid #ddd;
                                border-radius:12px;
                                padding:12px 16px;
                                color:#666;
                                font-size:14px;
                            ">
                                ⚪ No transcription factor impact data available for this LOF variant.
                            </div>
                            """,
                            unsafe_allow_html=True
                        )


    with tab2:
        # Load and combine all split files
        data_path = "./data/whole_genome_prediction_data/"  # change to the folder where your CSVs are
        all_files = sorted(glob.glob(os.path.join(data_path, "WGP_*.csv")))
        # Combine all parts into one DataFrame
        df_list = [pd.read_csv(f) for f in all_files] if all_files else []
        combined_wgp_df = pd.concat(df_list, ignore_index=True) if df_list else pd.DataFrame()
        
        st.markdown(
            """
            <h1 style="font-size:20px; font-weight:bold; color:#1f2937;">
                Enhancer regions predicted by DNABERT-Enhancer in human genome
            </h1>
            """,
            unsafe_allow_html=True
        )
        col1, col2 = st.columns([1, 4], gap="medium")
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
            if "filter_key_tab2" not in st.session_state:
                st.session_state.filter_key_tab2 = 0

            chrom_options = ["All"]
            if "chromosome" in combined_wgp_df.columns:
                chrom_options += sorted(combined_wgp_df["chromosome"].dropna().unique())
            selected_chrom = st.selectbox(
                "Chromosome",
                chrom_options,
                key=f"tab2_chrom_{st.session_state.filter_key_tab2}"
            )
            
            filtered_df = combined_wgp_df.copy()
            if selected_chrom != "All" and "chromosome" in filtered_df.columns:
                filtered_df = filtered_df[filtered_df["chromosome"] == selected_chrom]
            # Reset filters
            if st.button("🔄 Reset Filters", key="tab2_reset"):
                st.session_state.filter_key_tab2 += 1
                st.rerun()

        with col2:
            st.dataframe(filtered_df.head(50), use_container_width=True, height=500, hide_index=True)

            # Download option
            csv = filtered_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Download Filtered Enhancers (CSV)",
                data=csv,
                file_name="filtered_enhancer_regions.csv",
                mime="text/csv"
            )
        
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
            <li>Genome-wide enhancer predictions from DNABERT-Enhancer-350</li>
            <li>Candidate GOF/LOF variants predicted by DNABERT-Enhancer-350</li>
            <li>Interactive exploration by predicted effect, element class, chromosome, and more</li>
            <li>Transcription factor analysis</li>
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
    We thank all members of the Davuluri lab ( 548 The State University of the State University of New York at Stony Brook) and Dante Bolzan 
    (Ay Lab - La Jolla Institute for Immunology) for critical discussions and helpful 
    advice. This work was financially supported by grants from National Library of Medicine/National
    Institutes of Health funding – [R01LM01372201 to R.D., R35GM128938 to F.A]<br>
    For inquiries or collaborations, please contact: Tel: +1 (631) 638-2590; Email: Ramana.Davuluri@stonybrookmedicine.edu.<br>
    <br>
    Built with Streamlit 🎈
    </div>
        """, unsafe_allow_html=True)
